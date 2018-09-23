from peewee import *

database = SqliteDatabase('model-data.db', **{})

class UnknownField(object):
    def __init__(self, *_, **__): pass

class BaseModel(Model):
    class Meta:
        database = database

class CellType(BaseModel):
    acronym = TextField(column_name='Acronym', unique=True)
    id = AutoField(column_name='ID')
    name = TextField(column_name='Name', unique=True)

    class Meta:
        table_name = 'cell_type'

class CellLocation(BaseModel):
    cell_type = ForeignKeyField(column_name='Cell_Type_ID', field='id', model=CellType)
    id = AutoField(column_name='ID')
    x = FloatField()
    y = FloatField()

    class Meta:
        table_name = 'cell_location'

class GlomeruliLocation(BaseModel):
    enabled = IntegerField(column_name='Enabled')
    id = AutoField(column_name='ID')
    radius = FloatField(column_name='Radius')
    x = FloatField()
    y = FloatField()

    class Meta:
        table_name = 'glomeruli_location'

class Layer(BaseModel):
    acronym = TextField(column_name='Acronym', null=True)
    depth_order = IntegerField(column_name='Depth_Order')
    id = AutoField(column_name='ID')
    name = TextField(column_name='Name', unique=True)

    class Meta:
        table_name = 'layer'

class LayerVertex(BaseModel):
    id = AutoField(column_name='ID')
    layer = ForeignKeyField(column_name='Layer_ID', field='id', model=Layer)
    x = FloatField(null=True)
    y = FloatField(null=True)

    class Meta:
        table_name = 'layer_vertex'

class Parameter(BaseModel):
    name = TextField(primary_key=True)
    value = TextField()

    class Meta:
        table_name = 'parameter'

class SqliteSequence(BaseModel):
    name = UnknownField(null=True)  # 
    seq = UnknownField(null=True)  # 

    class Meta:
        table_name = 'sqlite_sequence'
        primary_key = False

