from sqlalchemy.sql.schema import Column
from sqlalchemy import Integer, String, Text, Date, Table
from database import Base

'''
class UnitInfo(Base):
    __tablename__ = "unit_info"

    unit_id = Column(String, primary_key=True)
    pdb_id = Column(String)
    chain = Column(String)
    chain_index = Column(Integer)

class LoopInfo(Base):
    __tablename__ = "loop_info"

    loop_id = Column(String, primary_key=True)
    unit_ids = Column(Text)
    loop_name = Column(Text)


class NrReleases(Base):
    __tablename__ = "nr_releases"

    nr_release_id = Column(String, primary_key=True)
    date = Column(Date)
    #classes = relationship("NrClasses", backref='nr_releases', lazy=True)


class NrClasses(Base):
    __tablename__ = "nr_classes"

    nr_class_id = Column(Integer, primary_key=True)
    name = Column(String)
    #nr_release_id = Column(String, ForeignKey("nr_releases.nr_release_id"))
    resolution = Column(String)
    #chains = relationship("NrChains", backref='nr_classes', lazy=True)


class NrChains(Base):
    __tablename__ = "nr_chains"

    nr_chain_id = Column(Integer, primary_key=True)
    ife_id = Column(String)
    nr_class_id = Column(Integer)
    #nr_class_id = Column(Integer, ForeignKey("nr_classes.nr_class_id"))
    #nr_release_id = Column(String, ForeignKey("nr_releases.nr_release_id"))


class UnitCorrespondence(Base):
    __tablename__ = "correspondence_units"

    unit_id_1 = Column(String, primary_key=True)
    unit_id_2 = Column(String)
    pdb_id_1 = Column(String)
    pdb_id_2 = Column(String)
    chain_name_2 = Column(String)


meta = Metadata()

test_view = Table("correspondence_units", meta,
                Column("unit_id_1", Integer, primary_key=True),
                autoload_with=engine
)  
'''  

# This is a view. Had to make each col as primary key for the query to work
class UnitCorrespondence(Base):
    __table__ = Table('correspondence_units', Base.metadata,
                Column('unit_id_1', String, primary_key=True),
                Column('chain_name_1', String, primary_key=True),
                Column('pdb_id_1', String, primary_key=True),
                Column('unit_id_2', String, primary_key=True),
                Column('chain_name_2', String, primary_key=True),
                Column('pdb_id_2', String, primary_key=True),
                )


class UnitRotation(Base):
    __table__ = Base.metadata.tables['unit_rotations']


class UnitCenter(Base):
    __table__ = Base.metadata.tables['unit_centers']


class NrChains(Base):
    __table__ = Base.metadata.tables['nr_chains']


class NrClasses(Base):
    __table__ = Base.metadata.tables['nr_classes']


class NrReleases(Base):
    __table__ = Base.metadata.tables['nr_releases']


class UnitInfo(Base):
    __table__ = Base.metadata.tables['unit_info']


class IfeInfo(Base):
    __table__ = Base.metadata.tables['ife_info']


class PDBInfo(Base):
    __table__ = Base.metadata.tables['pdb_info']


class LoopInfo(Base):
    __table__ = Base.metadata.tables['loop_info']


class UnitPairsInteractions(Base):
    __table__ = Base.metadata.tables['unit_pairs_interactions']


class UnitPairsDistances(Base):
    __table__ = Base.metadata.tables['unit_pairs_distances']


class ChainInfo(Base):
    __table__ = Base.metadata.tables['chain_info']
    
