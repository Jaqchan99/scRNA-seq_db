from sqlalchemy import create_engine, Column, String, DateTime, Text, Integer
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import sessionmaker
from datetime import datetime
import os

# SQLite 数据库配置
DATABASE_URL = "sqlite:///./scRNA_seq_platform.db"

engine = create_engine(DATABASE_URL, connect_args={"check_same_thread": False})
SessionLocal = sessionmaker(autocommit=False, autoflush=False, bind=engine)

Base = declarative_base()

class Submission(Base):
    """提交任务表"""
    __tablename__ = "submissions"
    
    id = Column(String, primary_key=True, index=True)
    status = Column(String, default="pending")  # pending, uploading, processing, completed, failed, validation_failed
    file_path = Column(String, nullable=True)
    report_path = Column(String, nullable=True)
    export_path = Column(String, nullable=True)
    error_message = Column(Text, nullable=True)
    created_at = Column(DateTime, default=datetime.now)
    updated_at = Column(DateTime, default=datetime.now, onupdate=datetime.now)

class MappingLog(Base):
    """映射日志表"""
    __tablename__ = "mapping_logs"
    
    id = Column(Integer, primary_key=True, index=True, autoincrement=True)
    submission_id = Column(String, index=True)
    step = Column(String)  # validation, gene_mapping, cell_type_mapping, export
    status = Column(String)  # success, warning, error
    details = Column(Text, nullable=True)
    created_at = Column(DateTime, default=datetime.now)

def get_db():
    """获取数据库会话"""
    db = SessionLocal()
    try:
        yield db
    finally:
        db.close()

def init_db():
    """初始化数据库表"""
    Base.metadata.create_all(bind=engine)
    print("📊 数据库初始化完成")

