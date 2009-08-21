"""SQLAlchemy definitions for the BioSQL database of biological items.

This provides a non-Seq-based interface to the BioSQL database through python.
This is useful if you have non-seq items to put into BioSQL, and should also
expose things not touched on with the Biopython BioSQL interface. Eventually
this would be a good target for merging.

http://www.biosql.org/wiki/Main_Page

Useful URLs for declarative style:
    https://www.bitbucket.org/stephane/model2/src/tip/transifex/model.py
    http://www.sqlalchemy.org/docs/05/sqlalchemy_ext_declarative.html
"""

def _initialize(Base):
    from sqlalchemy.orm import relation, mapper, dynamic_loader
    from sqlalchemy import MetaData, Table, Column, ForeignKey, Sequence
    from sqlalchemy import String, Unicode, Integer, DateTime, Float
    
    # -- Standard BioSQL tables
    
    class Biodatabase(Base):
        """Entry point to BioSQL databases.
        """
        __tablename__ = 'biodatabase'
        __table_args__ = {'mysql_engine':'InnoDB', 'autoload' : True}
        biodatabase_id = Column(Integer, primary_key = True)
        entries = relation("Bioentry", backref = "biodb")

    class Bioentry(Base):
        """The main bioentry object in BioSQL, containing a biological item.
        """
        __tablename__ = 'bioentry'
        __table_args__ = {'mysql_engine':'InnoDB', 'autoload' : True}
        bioentry_id = Column(Integer, primary_key = True)
        biodatabase_id = Column(Integer,
                ForeignKey('biodatabase.biodatabase_id'))
        qualifiers = relation("BioentryQualifierValue", backref = "bioentry")
        parent_maps = relation("BioentryRelationship", primaryjoin =
          "Bioentry.bioentry_id == BioentryRelationship.object_bioentry_id",
          lazy="dynamic")
        child_maps = relation("BioentryRelationship", primaryjoin =
          "Bioentry.bioentry_id == BioentryRelationship.subject_bioentry_id",
          order_by = "BioentryRelationship.object_bioentry_id.asc()",
          lazy="dynamic")
        features = relation("SeqFeature", backref="bioentry")
        sequence = relation("Biosequence")

    class Biosequence(Base):
        """Represent a sequence attached to a bioentry.
        """
        __tablename__ = 'biosequence'
        __table_args__ = {'mysql_engine':'InnoDB', 'autoload' : True}
        bioentry_id = Column(Integer, ForeignKey('bioentry.bioentry_id'),
                primary_key = True)
    
    class Ontology(Base):
        """Defined a high level dictionary of ontology key terms.
        """
        __tablename__ = 'ontology'
        __table_args__ = {'mysql_engine':'InnoDB', 'autoload' : True}
        ontology_id = Column(Integer, primary_key = True)
    
    class Term(Base):
        """Explicitly describe terminology used in key/value pair relationships
        """
        __tablename__ = 'term'
        __table_args__ = {'mysql_engine':'InnoDB', 'autoload' : True}
        term_id = Column(Integer, primary_key = True)
        ontology_id = Column(Integer, ForeignKey('ontology.ontology_id'))
        ontology = relation("Ontology", backref = "terms")

    class BioentryQualifierValue(Base):
        """A key/value annotation pair associated with a Bioentry.
        """
        __tablename__ = 'bioentry_qualifier_value'
        __table_args__ = {'mysql_engine':'InnoDB', 'autoload' : True}
        bioentry_id = Column(Integer,
                ForeignKey('bioentry.bioentry_id'), primary_key = True)
        term_id = Column(Integer, ForeignKey('term.term_id'), primary_key = True)
        rank = Column(Integer, primary_key = True)
        term = relation("Term", lazy=True)

    class BioentryRelationship(Base):
        """Define a relationship between two bioentry objects.
        """
        __tablename__ = 'bioentry_relationship'
        __table_args__ = {'mysql_engine':'InnoDB', 'autoload' : True}
        object_bioentry_id = Column(Integer,
                ForeignKey('bioentry.bioentry_id'), primary_key = True)
        subject_bioentry_id = Column(Integer, 
                ForeignKey('bioentry.bioentry_id'))
        term_id = Column(Integer, ForeignKey('term.term_id'),
                primary_key = True)
        rank = Column(Integer, primary_key = True)
        
        term = relation("Term")
        parent = relation("Bioentry", primaryjoin = 
          "Bioentry.bioentry_id == BioentryRelationship.object_bioentry_id")
        child = relation("Bioentry", primaryjoin =
          "Bioentry.bioentry_id == BioentryRelationship.subject_bioentry_id")

    seqfeature_dbxref_table = Table('seqfeature_dbxref', Base.metadata,
            Column('seqfeature_id', Integer, 
                ForeignKey('seqfeature.seqfeature_id')),
            Column('dbxref_id', Integer, ForeignKey('dbxref.dbxref_id')),
            Column('rank', Integer))
    
    class DBXref(Base):
        """Database cross reference.
        """
        __tablename__ = 'dbxref'
        __table_args__ = {'mysql_engine':'InnoDB', 'autoload' : True}
        dbxref_id = Column(Integer, primary_key = True)
    
    class SeqFeature(Base):
        """Provide a feature connected to a bioentry.
        """
        __tablename__ = 'seqfeature'
        __table_args__ = {'mysql_engine':'InnoDB', 'autoload' : True}
        seqfeature_id = Column(Integer, primary_key = True)
        bioentry_id = Column(Integer, ForeignKey('bioentry.bioentry_id'))
        type_term_id = Column(Integer, ForeignKey('term.term_id'))
        source_term_id = Column(Integer, ForeignKey('term.term_id'))
        type_term = relation("Term", primaryjoin =
            "SeqFeature.type_term_id == Term.term_id")
        source_term = relation("Term", primaryjoin =
            "SeqFeature.source_term_id == Term.term_id")
        qualifiers = relation("SeqFeatureQualifierValue")
        locations = relation("Location")
        dbxrefs = relation("DBXref", secondary=seqfeature_dbxref_table,
                order_by=seqfeature_dbxref_table.columns.rank)

    class SeqFeatureQualifierValue(Base):
        """A key/value annotation pair associated with a SeqFeature.
        """
        __tablename__ = 'seqfeature_qualifier_value'
        __table_args__ = {'mysql_engine':'InnoDB', 'autoload' : True}
        seqfeature_id = Column(Integer,
                ForeignKey('seqfeature.seqfeature_id'), primary_key = True)
        term_id = Column(Integer, ForeignKey('term.term_id'),
                primary_key = True)
        rank = Column(Integer, primary_key = True)
        term = relation("Term", lazy=True)
        
    class Location(Base):
        """Describe a location on a biological sequence.
        """
        __tablename__ = 'location'
        __table_args__ = {'mysql_engine':'InnoDB', 'autoload' : True}
        location_id = Column(Integer, primary_key = True)
        seqfeature_id = Column(Integer, ForeignKey('seqfeature.seqfeature_id'))
        dbxref_id = Column(Integer, ForeignKey('dbxref.dbxref_id'))
        qualifiers = relation("LocationQualifierValue")
        dbxref = relation("DBXref")
    
    class LocationQualifierValue(Base):
        """A key/value annotation pair associated with a Location.
        """
        __tablename__ = 'location_qualifier_value'
        __table_args__ = {'mysql_engine':'InnoDB', 'autoload' : True}
        location_id = Column(Integer,
                ForeignKey('location.location_id'), primary_key = True)
        term_id = Column(Integer, ForeignKey('term.term_id'),
                primary_key = True)
        rank = Column(Integer, primary_key = True)
        term = relation("Term", lazy=True)

    # ugly assignment of classes to the top level for use
    globals()['Biodatabase'] = Biodatabase
    globals()['Bioentry'] = Bioentry
    globals()['Biosequence'] = Biosequence
    globals()['BioentryQualifierValue'] = BioentryQualifierValue
    globals()['Ontology'] = Ontology
    globals()['Term'] = Term
    globals()['BioentryRelationship'] = BioentryRelationship
    globals()['SeqFeatureQualifierValue'] = SeqFeatureQualifierValue
    globals()['SeqFeature'] = SeqFeature
    globals()['DBXref'] = DBXref
    globals()['LocationQualifierValue'] = LocationQualifierValue
    globals()['Location'] = Location
