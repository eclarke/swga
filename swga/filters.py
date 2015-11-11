"""Provides filters for the `swga filter` step."""
from functools import wraps
import swga
import swga.primers

from swga.database import Primer

def _filter(fn):
    """Wrapper for filter methods.

    Updates the internal count of primers and ensures the statement
    executes after each filter method.
    """
    @wraps(fn)
    def func(self, *args, **kwargs):
        fn(self, *args, **kwargs)
        self._update_n()
        if not isinstance(self.primers, list):
            self.primers = list(self.primers)
        return self
    return func


def _update(fn):
    """Wrapper for update methods.

    Updates the database in chunks after each update method.
    """
    @wraps(fn)
    def func(self, *args, **kwargs):
        fn(self, *args, **kwargs)
        self.update_in_chunks
        return self
    return func


class Primers(object):
    """A list of primers with methods that operate on them."""

    def __init__(self, primers=None):
        """Create a new list of primers.

        :param primers: a list of Primer objects or a list of primer sequences.
        If None, selects all primers.
        """
        if primers is None:
            self.primers = Primer.select()
        else:
            self.primers = Primer.select().where(Primer.seq << primers)
    
        self.n = self.primers.count()

    @_filter
    def filter_min_fg_rate(self, min_bind):
        """
        Removes primers that bind less than the given rate to the foreground 
        genome.
        """
        self.primers = Primer.select().where(
            (Primer.seq << self.primers) &
            (Primer.fg_freq >= min_bind))

        swga.message(
            "{}/{} primers bind the foreground genome >= {} times"
            .format(self.primers.count(), self.n, min_bind))

    @_filter
    def filter_max_bg_rate(self, rate):
        """
        Removes primers that bind more than the given number of times to
        the background genome.
        """
        self.primers = Primer.select().where(
            (Primer.seq << self.primers) &
            (Primer.bg_freq <= rate))

        swga.message(
            "{}/{} primers bind the background genome <= {} times"
            .format(self.primers.count(), self.n, rate))


    @_filter
    def filter_tm_range(self, min_tm, max_tm):
        """
        Removes primers outside the given melting temp range.
        Finds melting temperatures if not already present.
        """
        swga.primers.update_Tms(self.primers)

        self.primers = Primer.select().where(
            (Primer.seq << self.primers) &
            (Primer.tm <= max_tm) &
            (Primer.tm >= min_tm))

        swga.message(
            "{}/{} primers have a melting temp between {} and {} C"
            .format(self.primers.count(), self.n, min_tm, max_tm))


    @_filter
    def limit_to(self, n):
        """
        Sorts by background binding rate, selects the top `n` least frequently
        binding, then returns those ordered by descending bg/fg ratio.
        """
        if n < 1:
            raise ValueError("n must be greater than 1")

        first_pass = (
            Primer.select().where(Primer.seq << self.primers)
            .order_by(Primer.bg_freq)
            .limit(n))

        self.primers = (
            Primer.select().where(Primer.seq << first_pass)
            .order_by(Primer.ratio.desc()))


    @_filter
    def filter_max_gini(self, gini_max, fg_genome_fp):
        """
        Filters out primers with Gini coefficients greater than the 
        max provided. Finds binding locations and Gini coefficients for primers 
        that do not have them already.

        :param gini_max: max Gini coefficient (0-1)
        """
        if 0 > gini_max > 1:
            raise ValueError("Gini coefficient must be between 0-1")

        (self
         .update_locations(fg_genome_fp)
         .update_gini(fg_genome_fp))

        self.primers = Primer.select().where(
            (Primer.seq << self.primers) &
            (Primer.gini <= gini_max))

        swga.message(
            "{}/{} primers have a Gini coefficient <= {}"
            .format(self.primers.count(), self.n, gini_max))


    def summarize(self):
        """Output the number of primers currently in list."""
        swga.message("{} primers satisfy all filters so far.".format(self.n))
        return self


    def activate(self, max_active=0):
        """Activates all the primers in the list.

        :param max_active: The maximum number expected to activate. Warns if
        fewer than this number.
        """
#        self.primers = self.primers.execute()
        n = (Primer.update(active=True)
            .where(Primer.seq << self.primers)
            .execute())

        swga.message("Marked {} primers as active.".format(n))
        if n < max_active:
            swga.warn(
                "Fewer than {} primers were selected ({} passed all the "
                "filters). You may want to try less restrictive filtering "
                "parameters.".format(max_active, n))
        return self


    @_update
    def update_locations(self, fg_genome_fp):
        """Finds binding locations for any primers that don't have them."""
        targets = list(Primer.select().where(
            (Primer.seq << self.primers) &
            (Primer._locations >> None)))

        for primer in targets:
            primer._update_locations(fg_genome_fp)


    @_update
    def update_gini(self, fg_genome_fp):
        """Calculates Gini coefs for any primers that don't have it."""
        targets = list(Primer.select().where(
            (Primer.seq << self.primers) &
            (Primer.gini >> None)))

        for primer in targets:
            primer._update_gini(fg_genome_fp)


    @_update
    def update_melt_temps(self):
        """Calculates melting temps for any primers that don't have it."""
        targets = list(Primer.select().where(
            (Primer.seq << self.primers) &
            (Primer.tm >> None)))

        for primer in targets:
            primer.update_tm()


    @_update
    def assign_ids(self):
        """Assigns sequential ids to active primers. 

        Resets any ids previously set.
        """
        # Reset all the primer IDs
        Primer.update(_id = -1).execute()
        
        primers = list(
            Primer.select()
            .where(Primer.seq << self.primers)
            .order_by(Primer.ratio.desc()).execute())

        for i, primer in enumerate(primers):
            primer._id = i + 1


    def _update_in_chunks(self, chunksize=100, show_progress=True, label="Updating primer db..."):
        """Updates the records for the list of primers in chunks.

        Technically, this performs an upsert, where primer records are updated 
        or inserted as needed. This is faster than updating the records for each
        primer individually.
        """

        def upsert_chunk(chunk):
            Primer.delete().where(Primer.seq << chunk).execute()
            Primer.insert_many(p.to_dict() for p in chunk).execute()

        swga.utils.chunk_iterator(
            itr=self.primers, 
            fn=upsert_chunk, 
            n=chunksize,
            show_progress=show_progress,
            label=label)


    def _update_n(self):
        n = self.primers.count()
        if n == 0:
            swga.error("No primers meet the given criteria.", exception=False)
        else:
            self.n = n
        

