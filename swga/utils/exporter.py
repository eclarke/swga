import csv
import peewee as pw
import collections
from playhouse.shortcuts import ManyToManyField

def export(model, rows, outfile, header=True):
    exported_fields = model.exported_fields()
    instance_fields = model.fields()
    writer = csv.DictWriter(
        outfile,
        exported_fields,
        delimiter="\t",
        extrasaction='ignore')
    
    if header:
        writer.writeheader()

    for row in rows:
        d = {}
        for field_name in row.exported_fields():
            field = getattr(row, field_name)
            if (not isinstance(field, basestring) and
                isinstance(field, collections.Iterable)):
                try:
                    field = ",".join([p.seq for p in field])
                except AttributeError:
                    field = ",".join([s._id for s in field])
            d[field_name] = field
        writer.writerow(d)
    
