import os
import swga
from swga import quote, __version__
import swga.database
from swga.commands._command import Command
from swga.database import Primer, Set
from swga.utils import fmtkv
import click
from peewee import fn

summary_template = '''
{header}
---------------
  Foreground genome:
    {fg_genome_fp}
  Background genome:
    {bg_genome_fp}
  Excluded sequences:
    {exclude_fp}


PRIMER SUMMARY
---------------
There are {nprimers} primers in the database. {if_no_primers_msg}

{nactive} are marked as active (i.e., they passed filter steps and will be
used to find sets of compatible primers.) {if_no_active_primers_msg}

The average number of foreground genome binding sites is {avg_fg_bind:.0f}.
    (avg binding / genome_length = {fg_bind_ratio:05f})
The average number of background genome binding sites is {avg_bg_bind:.0f}.
    (avg binding / genome_length = {bg_bind_ratio:05f})

{melting_tmp_msg}


SETS SUMMARY
---------------
There are {nsets} sets in the database.
{set_msg}\
'''

best_set_desc = '''\
The best scoring set is #{best_set}, with {bs_size} primers and a score of {bs_score:03f}.
Various statistics:
 {bs_stats}
The primers in Set {best_set} are:
{bs_primers}
'''


class Summary(Command):

    def run(self):
        self.summary_msg = summary_template
        self.best_set_desc = best_set_desc
        avg_fg_bind, avg_bg_bind, nprimers = (
            Primer
            .select(fn.Avg(Primer.fg_freq),
                    fn.Avg(Primer.bg_freq),
                    fn.Count(Primer.seq))
            .scalar(as_tuple=True))

        if (avg_fg_bind is None) or (avg_bg_bind is None):
            (avg_fg_bind, avg_bg_bind) = (0, 0)

        fg_bind_ratio = avg_fg_bind / float(self.fg_length)
        bg_bind_ratio = avg_bg_bind / float(self.bg_length)
        nactive = Primer.select().where(Primer.active == True).count()

        min_tm, max_tm, avg_tm = (
            Primer
            .select(fn.Min(Primer.tm),
                    fn.Max(Primer.tm),
                    fn.Avg(Primer.tm))
            .where(Primer.active == True)
            .scalar(as_tuple=True))

        nsets = Set.select(fn.Count(Set._id)).scalar()

        if nsets > 0:
            bs = Set.select().order_by(Set.score).limit(1).get()
            bs_primers = ", ".join(
                swga.database.get_primers_for_set(bs._id)).strip()
            best_set = bs._id
            bs_size = bs.set_size
            bs_score = bs.score
            bs_stats = "- " + "\n - ".join(
                fmtkv(k, bs.__dict__['_data'][k])
                for k in bs.exported_fields()
                if k not in ["_id", "pids", "score", "primers"]
            )
            self.best_set_desc = self.best_set_desc.format(**locals())

        if_no_primers_msg = click.style(
            "Run `swga count` to find possible primers."
            if nprimers == 0 else "", fg='green')
        if_no_active_primers_msg = click.style(
            "Run `swga filter` to identify primers to use."
            if nactive == 0 else "", fg='green')
        melting_tmp_msg = (
            "The melting temp of the primers ranges between {min_tm:.2f}C and "
            "{max_tm:.2f}C with an average of {avg_tm:.2f}C."
            if nactive > 0 and min_tm and max_tm else
            "No melting temps have been calculated yet.").format(**locals())
        ifzero_sets_msg = click.style(
            "Run `swga find_sets` after identifying valid primers to begin "
            "collecting sets.\n", fg='green')

        set_msg = (self.best_set_desc if nsets > 0 else ifzero_sets_msg)

        primer_db = os.path.abspath(self.primer_db)
        nprimers = click.style(str(nprimers), bold=True, fg='blue')
        nactive = click.style(str(nactive), bold=True, fg='blue')
        nsets = click.style(str(nsets), bold=True, fg='blue')

        self.header = click.style("swga v{}".format(__version__), fg='green')

        # Copy all the relevant values into one dict
        values = self.__dict__.copy()
        values.update(locals())

        # Format the summary message with all the calculated values

        self.summary_msg = self.summary_msg.format(**values)
        click.echo(quote(self.summary_msg, quote="  ", width=200))
