import cobra
import pandas as pd
import pymysql
from sqlalchemy import create_engine

with open('sql_credentials', 'r') as fid:
    sql_connection_string = fid.read()

engine = create_engine(sql_connection_string, echo=False)

sql_query = """
SELECT 
    t2.eid1,
    t2.iid1 as enzyme_id,
    t1.strV,
    t3.eid1,
    t3.iid1 as publication_id,
    t5.name,
    t3.row1,
    t3.eid3,
    t3.value_type,
    t3.intV,
    t3.floatV,
    t3.strV
FROM
    entity_instance_entity t0,
    `pw27-27`.instance_instance_instance t2,
    instance_property t1,
    instance_property t3,
    property t5
WHERE
    t0.eid1 = 13 AND t0.eid2 = 16
        AND t0.iid2 = 2
        AND t1.eid1 = t0.eid3
        AND t1.pid = 1
        AND t1.eid3 = 3
        AND t1.strv LIKE '%%{ec_number})%%'
        AND t2.eid1 = t1.eid1
        AND t2.iid1 = t1.iid1
        AND t2.eid2 = 17
        AND t2.iid2 = 1
        AND t2.eid3 = 35
        AND t3.eid1 = t2.eid3
        AND t3.iid1 = t2.iid3
        AND t5.pid = t3.pid
        AND t5.name = 'kcatValue'

"""

# df = pd.read_sql_query(sql_query.format(ec_number='1.2.4.4'), con = engine)
# grouped = df.groupby('publication_id')

ecoli_with_ec = cobra.io.load_json_model('../models/iJO1366_with_xrefs.json')

reaction_ecs = dict()

EC_FIELD = 'ec_numbers'

for r in ecoli_with_ec.reactions:
    if EC_FIELD in r.notes:
        ec_numbers = r.notes[EC_FIELD]

        for ec in ec_numbers:
            if ec not in reaction_ecs:
                df = pd.read_sql_query(sql_query.format(ec_number=ec), con = engine)
                the_value = df[df['value_type'] == 2]['floatV'].max()
                reaction_ecs[ec]  = the_value

ec_df = pd.DataFrame.from_dict(reaction_ecs, orient='index')
ec_df.columns = ['kcat']
ec_df = ec_df[ec_df > 0]
ec_df.dropna().to_csv('info_ecoli/aggregated_kcats.csv')