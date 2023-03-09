import psycopg2
import argparse
import psycopg2
import psycopg2.extras


def get_args():
    args = argparse.ArgumentParser()

    args.add_argument("server")
    args.add_argument("database")
    args.add_argument("user")
    args.add_argument("password")

    args.add_argument("contig_name")

    return args

def main():
    args = vars(get_args().parse_args())


    conn_string = "host='%(server)s' dbname='%(database)s' user='%(user)s' password='%(password)s'" % args

    conn = psycopg2.connect(conn_string)

    cursor = conn.cursor(cursor_factory=psycopg2.extras.DictCursor)

    print "Dropping strains_meta...",
    cursor.execute("DROP TABLE IF EXISTS strain_meta")
    print "Done"

    print "Adding ignored_pos...",
    cursor.execute("CREATE TABLE ignored_pos (id integer PRIMARY KEY, contig character varying, pos integer)")
    print "Done"

    print "Dropping locus_tag from variants...",
    try:
        cursor.execute("ALTER TABLE variants DROP COLUMN locus_tag")
        print "Done"
    except:
        conn.rollback()
        cursor.execute("DROP TABLE IF EXISTS strain_meta")
        cursor.execute("CREATE TABLE ignored_pos (id integer PRIMARY KEY, contig character varying, pos integer)")
        print "Skip"

    print "Adding contig to variants...",
    cursor.execute("ALTER TABLE variants ADD COLUMN contig varchar")
    cursor.execute("UPDATE variants SET contig=%(contig_name)s", args)

    print "Done"

    print "Dropping time_of_upload from strains_snps...",
    cursor.execute("ALTER TABLE strains_snps DROP COLUMN time_of_upload")

    print "Done"

    print "Filling ignored_pos table ...",
    cursor.execute("SELECT ignored_pos FROM strains_snps")
    data = set()
    for row in cursor:
        data.update([ (pos, args["contig_name"], pos) for pos in row["ignored_pos"] ] )
    cursor.executemany("INSERT INTO ignored_pos (id, contig, pos) VALUES (%s,%s,%s)", data)
    print "Done"

    cursor.execute("CREATE SEQUENCE ignored_pos_id_seq")
    cursor.execute("SELECT setval('ignored_pos_id_seq', (select max(id) from ignored_pos));")
    cursor.execute("ALTER TABLE ignored_pos ALTER COLUMN id SET DEFAULT nextval('ignored_pos_id_seq');")

    conn.commit()
    conn.close()

    return 0

if __name__ == "__main__":
    exit(main())
