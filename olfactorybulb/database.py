from peewee import *

database = SqliteDatabase('olfactorybulb/model-data.sqlite', **{})

class UnknownField(object):
    def __init__(self, *_, **__): pass

class BaseModel(Model):
    class Meta:
        database = database

class CellType(BaseModel):
    id = TextField(column_name='ID', primary_key=True)
    name = TextField(column_name='Name', unique=True)

    class Meta:
        table_name = 'cell_type'

class Cell(BaseModel):
    id = AutoField(column_name='ID')
    type = ForeignKeyField(column_name='Type', field='id', model=CellType)
    x = FloatField()
    y = FloatField()
    z = FloatField()

    class Meta:
        table_name = 'cell'

class Glomerulus(BaseModel):
    id = AutoField(column_name='ID')
    radius = FloatField(column_name='Radius')
    x = FloatField()
    y = FloatField()
    z = FloatField()

    class Meta:
        table_name = 'glomerulus'

class Layer(BaseModel):
    depth_order = IntegerField(column_name='Depth_Order')
    id = TextField(column_name='ID', primary_key=True)
    name = TextField(column_name='Name', unique=True)

    class Meta:
        table_name = 'layer'

class LayerVertex(BaseModel):
    id = AutoField(column_name='ID')
    layer = ForeignKeyField(column_name='Layer_ID', field='id', model=Layer)
    x = FloatField()
    y = FloatField()
    z = FloatField()

    class Meta:
        table_name = 'layer_vertex'

class LfpElectrode(BaseModel):
    id = AutoField(column_name='ID')
    label = TextField(column_name='Label')
    x = FloatField()
    y = FloatField()
    z = FloatField()

    class Meta:
        table_name = 'lfp_electrode'

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

