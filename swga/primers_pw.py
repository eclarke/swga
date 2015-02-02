import peewee as pw

db = pw.SqliteDatabase(None)

class Primer(pw.Model):
    seq = pw.TextField(unique=True)
    fg_freq = pw.IntegerField(default=0)
    bg_freq = pw.IntegerField(default=0)
    ratio = pw.FloatField(default=0.0)
    tm = pw.FloatField(default=0.0)
    locations = pw.TextField(default="")

    class Meta:
        database = db


    
    
