"""Provides filters for the `swga filter` step."""
import swga
import swga.primers

from swga.database import Primer

class Filter():

    def __init__(self, primers=None, seqs=None):
        """Build a new filter based on the provided primers."""
        if seqs:
            self.primers = Primer.select().where(Primer.seq << seqs)
        elif primers:
            self.primers = primers
        else:
            raise ValueError("Must provide either Primers or primer sequences")


    def summarize(self):
        swga.message(
            "{} primers satisfy all filters so far."
            .format(self.primers.count()))
        return self


    def min_fg_rate(self, rate):
        """
        Removes primers that bind less than the given rate to the foreground 
        genome.
        """
        self.primers = Primer.select().where(
            (Primer.seq << self.primers) &
            (Primer.fg_freq >= rate))
        swga.message(
            "{} primers bind the foreground genome >= {} times"
            .format(self.primers.count(), rate))
        return self


    def max_bg_rate(self, rate):
        """
        Removes primers that bind more than the given number of times to
        the background genome.
        """
        self.primers = Primer.select().where(
            (Primer.seq << self.primers) &
            (Primer.bg_freq <= rate))
        swga.message(
            "{} primers bind the background genome <= {} times"
            .format(self.primers.count(), rate))
        return self


    def tm_range(self, min_tm, max_tm):
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
            "{} primers have a melting temp between {} and {} C"
            .format(self.primers.count(), min_tm, max_tm))
        return self


    def limit_ratio(self, n):
        """
        Sorts by background binding rate, selects the top `n` least frequently
        binding, then returns those ordered by descending bg/fg ratio.
        """
        first_pass = (
            Primer.select().where(Primer.seq << self.primers)
            .order_by(Primer.bg_freq)
            .limit(n)
        )
        self.primers = (
            Primer.select().where(Primer.seq << first_pass)
            .order_by(Primer.ratio.desc())
        )
        return self

    def max_gini(self, gini_max, fg_genome_fp):
        """
        Filters out primers with Gini coefficients greater than the 
        max provided. 
        Finds binding locations and Gini coefficients for primers that
        do not have them already.
        """
        swga.primers.update_locations(self.primers, fg_genome_fp)
        swga.primers.update_gini(self.primers, fg_genome_fp)
        self.primers = Primer.select().where(
            (Primer.seq << self.primers) &
            (Primer.gini <= gini_max))
        swga.message(
            "{} primers have a Gini coefficient <= {}"
            .format(self.primers.count(), gini_max))
        return self
        

