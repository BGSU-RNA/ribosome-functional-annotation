from sqlalchemy import create_engine
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import sessionmaker, scoped_session
from contextlib import contextmanager

'''
# you would need to change the username & password
SQLALCHEMY_DATABASE_URL = "mysql://sria:sria@127.0.0.1/rna3dhub-next"
engine = create_engine(SQLALCHEMY_DATABASE_URL)

session = sessionmaker(autocommit=False, autoflush=False, bind=engine)

Base = declarative_base()
'''

engine = create_engine('mysql://sria:sria@127.0.0.1/rna3dhub-next', convert_unicode=True, echo=False)
Base = declarative_base()
Base.metadata.reflect(engine)

SessionLocal = scoped_session(sessionmaker(bind=engine))

@contextmanager
def db_session():
	db = SessionLocal()
	try:
		yield db
		db.commit()
	except:
		db.rollback()
		raise
	finally:
		db.close()
