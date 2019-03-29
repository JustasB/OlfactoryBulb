from peewee import *

database = SqliteDatabase('olfactorybulb/model-data.sqlite', **{})

class UnknownField(object):
    def __init__(self, *_, **__): pass

class BaseModel(Model):
    class Meta:
        database = database

class CellType(BaseModel):
    id = TextField(primary_key=True)
    name = TextField(unique=True)

    class Meta:
        table_name = 'cell_type'

class Cell(BaseModel):
    type = ForeignKeyField(column_name='type', model=CellType)
    x = FloatField()
    y = FloatField()
    z = FloatField()

    class Meta:
        table_name = 'cell'

class Glomerulus(BaseModel):
    radius = FloatField()
    x = FloatField()
    y = FloatField()
    z = FloatField()

    class Meta:
        table_name = 'glomerulus'

class Layer(BaseModel):
    depth_order = IntegerField()
    id = TextField(primary_key=True)
    name = TextField(unique=True)

    class Meta:
        table_name = 'layer'

class LayerVertex(BaseModel):
    layer = ForeignKeyField(column_name='layer_id', model=Layer)
    x = FloatField()
    y = FloatField()
    z = FloatField()

    class Meta:
        table_name = 'layer_vertex'

class LfpElectrode(BaseModel):
    label = TextField()
    x = FloatField()
    y = FloatField()
    z = FloatField()

    class Meta:
        table_name = 'lfp_electrode'

class Source(BaseModel):
    journal_book_title = TextField(null=True)
    short_title = TextField()
    title = TextField(null=True)
    url = TextField()
    year = IntegerField(null=True)

    class Meta:
        table_name = 'source'

class Property(BaseModel):
    id = TextField(primary_key=True)
    is_emergent = IntegerField(null=True)
    mean = FloatField(null=True)
    n = IntegerField(null=True)
    name = TextField(null=True)
    type = TextField(null=True)
    picked_value = FloatField(null=True)
    std = FloatField(null=True)

    class Meta:
        table_name = 'property'

class Measurement(BaseModel):
    mean = FloatField(null=True)
    n = IntegerField()
    notes = TextField(null=True)
    property = ForeignKeyField(column_name='property_id', field='id', model=Property)
    source = ForeignKeyField(column_name='source_id', field='id', model=Source)
    std = FloatField(null=True)

    class Meta:
        table_name = 'measurement'

class SqliteSequence(BaseModel):
    name = UnknownField(null=True)  # 
    seq = UnknownField(null=True)  # 

    class Meta:
        table_name = 'sqlite_sequence'
        primary_key = False

