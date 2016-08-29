##########################################################################################################
#Class for forms used on input page
##########################################################################################################
from wtforms import Form, FileField, StringField, IntegerField, BooleanField
from wtforms.validators import DataRequired

class RnaForm(Form):
    #required input files
    ref_fasta = FileField('ref_fasta')
    
    #optional input files
    ref_shape = FileField('ref_shape')
    ref_db = FileField('ref_db')
    map_fasta = FileField('map_fasta')
    map_db = FileField('map_db')

    #advanced input
    ref_size = IntegerField('ref_size', default=10)
    ref_rgstart = IntegerField('ref_rgstart', default=1)
    ref_rgend = IntegerField('ref_rgend', default=None)
    ref_maxd = IntegerField('ref_maxd', default=None)
    email = StringField('email', default=None)
