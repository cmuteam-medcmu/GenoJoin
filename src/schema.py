from sqlalchemy import (
    BigInteger,
    Column,
    ForeignKey,
    Integer,
    String,
    Text,
    UniqueConstraint,
)
from sqlalchemy.orm import declarative_base, relationship

# Base class
Base = declarative_base()


# Variant table
class Variant(Base):
    __tablename__ = "variants"
    id = Column(Integer, primary_key=True, autoincrement=True)
    chrom = Column(Text, nullable=False)
    pos = Column(BigInteger, nullable=False)
    ref = Column(Text, nullable=False)
    alt = Column(Text, nullable=False)
    samples = relationship("Sample", back_populates="variant")

    __table_args__ = (
        UniqueConstraint("chrom", "pos", "ref", "alt", name="uix_variant_composite"),
    )


# Sample table
class Sample(Base):
    __tablename__ = "samples"
    id = Column(Integer, primary_key=True, autoincrement=True)
    name = Column(String, nullable=False)
    vid = Column(Integer, ForeignKey("variants.id"), nullable=False)
    pat = Column(String, nullable=False)
    variant = relationship("Variant", back_populates="samples")
