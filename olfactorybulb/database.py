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

class Source(BaseModel):
    error_type = TextField(null=True)
    id = TextField(primary_key=True)
    journal_book_title = TextField(null=True)
    liquid_junction_potential = TextField(null=True)
    notes = TextField(null=True)
    publication_class = TextField(null=True, unique=True)
    temperature_celsius = UnknownField(null=True)  # FLOAT
    title = TextField(null=True)
    url = TextField()

    class Meta:
        table_name = 'source'

class CellModel(BaseModel):
    rowid = IntegerField(primary_key=True)
    cell_type = ForeignKeyField(column_name='cell_type', field='id', model=CellType, null=True)
    isolated_model_class = TextField(null=True)
    name = TextField(null=True)
    source = ForeignKeyField(column_name='source_id', field='id', model=Source, null=True)

    class_name = TextField(null=True)
    apical_dendrite_start = IntegerField(null=True)
    apical_dendrite_end = IntegerField(null=True)
    apical_dendrite_reach = FloatField(null=True)
    tufted_dend_root = TextField(null=True)

    class Meta:
        table_name = 'cell_model'
        primary_key = False

class Glomerulus(BaseModel):
    rowid = IntegerField(primary_key=True)
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

class Property(BaseModel):
    id = TextField(primary_key=True)
    is_emergent = IntegerField(null=True)
    units = TextField(null=True)
    mean = FloatField(null=True)
    n = IntegerField(null=True)
    name = TextField(null=True)
    picked_value = FloatField(null=True)
    std = FloatField(null=True)
    test_class_generic = TextField(null=True)
    type = TextField(null=True)

    class Meta:
        table_name = 'property'

class Measurement(BaseModel):
    mean = FloatField(null=True)
    n = IntegerField()
    notes = TextField(null=True)
    property = ForeignKeyField(column_name='property_id', field='id', model=Property)
    source = ForeignKeyField(column_name='source_id', field='id', model=Source, null=True)
    std = FloatField(null=True)

    class Meta:
        table_name = 'measurement'
        primary_key = False

class Odor(BaseModel):
    rowid = IntegerField(primary_key=True)
    name = TextField(unique=True)

    class Meta:
        table_name = 'odor'

class OdorGlom(BaseModel):
    rowid = IntegerField(primary_key=True)
    glom = ForeignKeyField(column_name='glom_id', field='rowid', model=Glomerulus)
    odor = ForeignKeyField(column_name='odor_id', field='rowid', model=Odor)

    intensity = FloatField()

    class Meta:
        table_name = 'odor_glom'

class SqliteSequence(BaseModel):
    name = UnknownField(null=True)  # 
    seq = UnknownField(null=True)  # 

    class Meta:
        table_name = 'sqlite_sequence'
        primary_key = False

