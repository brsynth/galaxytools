from sqlalchemy import create_engine, text

db_uri = "postgresql://postgres:RK17@localhost:5432/test_fragments_db" #adapt with your URI's DB
engine = create_engine(db_uri)

with engine.connect() as conn:
    result = conn.execute(text("""
        SELECT fragment, sequence, annotation
        FROM sample
        ORDER BY fragment
    """))

    print("Full contents of fragments in DB:\n")
    for row in result:
        print(f" Fragment: {row.fragment}")
        print(" Sequence:")
        print(row.sequence)
        print("\n Annotation:")
        print(row.annotation)
        print("-" * 80)