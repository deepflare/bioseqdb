CREATE DOMAIN nucl_seq AS TEXT
CHECK(
   VALUE ~ '^[ACGTX]*$'
);

CREATE DOMAIN aa_seq AS TEXT
CHECK(
   VALUE ~ '^[ARNDCQEGHILKMFPSTWYV\*]*$'
);

CREATE DOMAIN amb_aa_seq AS TEXT
CHECK(
   VALUE ~ '^[ARNDCQEGHILKMFPSTWYVBJZX\*]*$'
);

-- TABLES STARTING WITH A PREFIX
CREATE OR REPLACE FUNCTION index_tables_for(
	preffix_t_name information_schema.sql_identifier
) RETURNS TABLE (
	"schema" information_schema.sql_identifier, 
	"tbl_name" information_schema.sql_identifier
) AS
$BODY$
BEGIN
	RETURN QUERY (
        SELECT
            "table_schema", "table_name"
        FROM
            information_schema.tables
        WHERE
            table_type = 'BASE TABLE'
        -- AND
        --    table_schema = changed_object.schema_name
        AND
            table_name LIKE (preffix_t_name || '%' || '__index')
	);
END;
$BODY$
LANGUAGE plpgsql;


CREATE OR REPLACE FUNCTION show_domain_usage(
	d_name text
) RETURNS SETOF information_schema.column_domain_usage AS
$BODY$
BEGIN
	RETURN QUERY (
		SELECT 
			*
		FROM 
			information_schema.column_domain_usage
		WHERE
			"domain_name" = d_name
	);
END;
$BODY$
LANGUAGE plpgsql;

CREATE OR REPLACE FUNCTION create_index_table(
	t_name information_schema.sql_identifier, 
	c_name information_schema.sql_identifier
) RETURNS void AS
$BODY$
BEGIN
	-- There could be a foreign key, but if someone wanted to delete 
	-- their table they would have to use 'CASCADE'
	-- FOREIGN KEY (seq_id) REFERENCES %I(id) 
	EXECUTE format(
		-- using text, so aa_seq, amb_aa_seq and nucl_seq can all be used
	    'CREATE TABLE %I (
	        kmer text NOT NULL,
	        starting_position INT NOT NULL, 
	        seq_id BIGINT NOT NULL
		)',
	   	(t_name || '_' || c_name || '__index'), 
	   	(t_name)
	);
	EXECUTE format(
	    'CREATE INDEX %I ON %I (kmer);',
	   	(t_name || '_' || c_name || '__index' || '_kmer_idx'), 
	   	(t_name || '_' || c_name || '__index')
	);
	EXECUTE format(
	    'CREATE INDEX %I ON %I (seq_id);',
	   	(t_name || '_' || c_name || '__index' || '_seq_id_idx'), 
	   	(t_name || '_' || c_name || '__index')
	);
END;
$BODY$
LANGUAGE plpgsql;

-- INSERT OCCURRANCE
CREATE OR REPLACE FUNCTION insert_index(
	t_name information_schema.sql_identifier, 
	c_name information_schema.sql_identifier, 

	kmer text,
	starting_position INT, 
	seq_id BIGINT
) RETURNS void AS
$BODY$
BEGIN 
	EXECUTE FORMAT('
		INSERT INTO %I("kmer", "starting_position", "seq_id") VALUES(%L, %L, %L)', 
	   (t_name || '_' || c_name || '__index'),
	   kmer, starting_position, seq_id
	);
END;
$BODY$
LANGUAGE plpgsql;

CREATE OR REPLACE FUNCTION insert_indexes(
	t_name information_schema.sql_identifier, 
	c_name information_schema.sql_identifier, 
	
	seq_id BIGINT,
	seq text 
) RETURNS void AS
$BODY$
DECLARE
	ch text;
	kmer text := '';
	ending_position SMALLINT := -1;

BEGIN
    FOR ch IN SELECT regexp_split_to_table(seq, '') LOOP
		ending_position := (ending_position + 1);
		kmer := (kmer || ch);
		
	 	IF ending_position >= 7 THEN 
			kmer := substr(kmer, 2); 
	  	END IF;
	
	  	IF ending_position >= 6 THEN 	  		
			EXECUTE insert_index(t_name, c_name, kmer, ending_position - 6, seq_id);
	  	END IF;

    END LOOP;
    RETURN;
END;
$BODY$
LANGUAGE plpgsql;

CREATE OR REPLACE FUNCTION update_index_trigger() RETURNS TRIGGER AS
$BODY$
DECLARE 
	c_name information_schema.sql_identifier;

	seq text;
BEGIN
	FOR c_name IN (
		SELECT column_name FROM show_domain_usage('nucl_seq')
		WHERE table_name = TG_TABLE_NAME 
		UNION
		SELECT column_name FROM show_domain_usage('aa_seq')
		WHERE table_name = TG_TABLE_NAME 
		UNION
		SELECT column_name FROM show_domain_usage('amb_aa_seq')
		WHERE table_name = TG_TABLE_NAME 
	)
	LOOP 
		EXECUTE FORMAT('DELETE FROM %I WHERE "seq_id" = %L', 
			(TG_TABLE_NAME || '_' || c_name || '__index'), OLD.id);

		EXECUTE FORMAT ('SELECT %I FROM %I WHERE id = %L', 
			c_name, TG_TABLE_NAME, NEW.id) INTO seq;

		EXECUTE FORMAT ('SELECT insert_indexes(%L, %L, %s, %L)', 
			TG_TABLE_NAME, c_name, NEW.id, seq);
		
	END LOOP;

	RETURN NEW;
END;
$BODY$
LANGUAGE plpgsql;

CREATE OR REPLACE FUNCTION insert_index_trigger() RETURNS TRIGGER AS
$BODY$
DECLARE 
	c_name information_schema.sql_identifier;

	seq text;
BEGIN
	FOR c_name IN (
		SELECT column_name FROM show_domain_usage('nucl_seq')
		WHERE table_name = TG_TABLE_NAME 
		UNION
		SELECT column_name FROM show_domain_usage('aa_seq')
		WHERE table_name = TG_TABLE_NAME 
		UNION
		SELECT column_name FROM show_domain_usage('amb_aa_seq')
		WHERE table_name = TG_TABLE_NAME 
	)
	LOOP 
		EXECUTE FORMAT ('SELECT %I FROM %I WHERE id = %L', 
			c_name, TG_TABLE_NAME, NEW.id) INTO seq;

		EXECUTE FORMAT ('SELECT insert_indexes(%L, %L, %s, %L)', 
			TG_TABLE_NAME, c_name, NEW.id, seq);
		
	END LOOP;

	RETURN NEW;
END;
$BODY$
LANGUAGE plpgsql;

CREATE OR REPLACE FUNCTION delete_index_trigger() RETURNS TRIGGER AS
$BODY$
DECLARE 
	c_name information_schema.sql_identifier;
BEGIN
	FOR c_name IN (
		SELECT column_name FROM show_domain_usage('nucl_seq')
		WHERE table_name = TG_TABLE_NAME 
		UNION
		SELECT column_name FROM show_domain_usage('aa_seq')
		WHERE table_name = TG_TABLE_NAME 
		UNION
		SELECT column_name FROM show_domain_usage('amb_aa_seq')
		WHERE table_name = TG_TABLE_NAME 
	)
	LOOP 
		
		EXECUTE FORMAT('
			DELETE FROM %I WHERE "seq_id" = %L
		', 
		(TG_TABLE_NAME || '_' || c_name || '__index'), 
		OLD.id);

	END LOOP;

	RETURN NEW;
END;
$BODY$
LANGUAGE plpgsql;

-- CREATE create_mmseqs_index_tables_after_change
CREATE OR REPLACE FUNCTION create_mmseqs_index_tables_after_change()
        RETURNS event_trigger LANGUAGE plpgsql AS $$
DECLARE
    row record;
    changed_object record;

	t_name information_schema.sql_identifier; 
	c_name information_schema.sql_identifier;

	seq_id BIGINT;
	seq text;
BEGIN
    FOR changed_object IN SELECT * FROM pg_event_trigger_ddl_commands()
    LOOP		
		-- Creating new tables with indexes for given table
    	IF changed_object.object_type = 'table' 
		AND changed_object.object_identity NOT LIKE '%__index' 
		THEN
			SELECT SPLIT_PART(changed_object.object_identity, '.', 2) INTO t_name;
			
			-- Deleting ALL previous tables
		    FOR row IN SELECT * FROM index_tables_for(SPLIT_PART(changed_object.object_identity, '.', 2))
		    LOOP EXECUTE format('DROP TABLE %I CASCADE', row.tbl_name); END LOOP;
			
			-- Creating new tables
			FOR c_name IN (
				SELECT column_name FROM show_domain_usage('nucl_seq')
				WHERE "table_schema" || '.' || "table_name" = changed_object.object_identity
				UNION
				SELECT column_name FROM show_domain_usage('aa_seq')
				WHERE "table_schema" || '.' || "table_name" = changed_object.object_identity
				UNION
				SELECT column_name FROM show_domain_usage('amb_aa_seq')
				WHERE "table_schema" || '.' || "table_name" = changed_object.object_identity
			)
			LOOP 

			   	EXECUTE create_index_table(SPLIT_PART(changed_object.object_identity, '.', 2), c_name);
												
				FOR seq_id, seq IN EXECUTE format('SELECT id, %I as clmn from %I', c_name, t_name) LOOP
					EXECUTE insert_indexes(t_name, c_name, seq_id, seq);
				END LOOP;
				
			END LOOP;

			EXECUTE FORMAT('
				DROP TRIGGER IF EXISTS %I
				ON %I;

				CREATE TRIGGER %I
			    AFTER UPDATE ON %I
			    FOR EACH ROW
			    WHEN (OLD IS DISTINCT FROM NEW)
			    EXECUTE PROCEDURE update_index_trigger()
			', 
			t_name || '__update_trigger', 
			t_name,
			t_name || '__update_trigger', 
			t_name
			);

			EXECUTE FORMAT('
				DROP TRIGGER IF EXISTS %I
				ON %I;

				CREATE TRIGGER %I
			    AFTER INSERT ON %I
			    FOR EACH ROW
			    EXECUTE PROCEDURE insert_index_trigger()
			', 
			t_name || '__insert_trigger', 
			t_name,
			t_name || '__insert_trigger', 
			t_name
			);

			EXECUTE FORMAT('
				DROP TRIGGER IF EXISTS %I
				ON %I;

				CREATE TRIGGER %I
			    AFTER DELETE ON %I
			    FOR EACH ROW
			    EXECUTE PROCEDURE delete_index_trigger()
			', 
			t_name || '__delete_trigger', 
			t_name,
			t_name || '__delete_trigger', 
			t_name
			);

		END IF;
    END LOOP;
END;
$$;

					
CREATE EVENT TRIGGER 
	indexise_after_change_trigger 
ON 
	ddl_command_end
EXECUTE PROCEDURE 
	create_mmseqs_index_tables_after_change();

CREATE OR REPLACE FUNCTION test_event_trigger_for_drops()
        RETURNS event_trigger LANGUAGE plpgsql AS $$
DECLARE
	row record;
    changed_object record;

	t_name information_schema.sql_identifier; 
BEGIN
    FOR changed_object IN SELECT * FROM pg_event_trigger_dropped_objects()
    LOOP
    	IF changed_object.object_type = 'table' AND changed_object.object_name NOT LIKE '%__index' THEN	                   
			SELECT SPLIT_PART(changed_object.object_identity, '.', 2) INTO t_name;
  
		    FOR row IN SELECT * FROM index_tables_for(SPLIT_PART(changed_object.object_identity, '.', 2))
		    LOOP EXECUTE format('DROP TABLE %I CASCADE', row.tbl_name); END LOOP;
	           
			EXECUTE FORMAT('DROP TRIGGER IF EXISTS %I ON %I', t_name || '__update_trigger', t_name);
			EXECUTE FORMAT('DROP TRIGGER IF EXISTS %I ON %I', t_name || '__insert_trigger', t_name);
			EXECUTE FORMAT('DROP TRIGGER IF EXISTS %I ON %I', t_name || '__delete_trigger', t_name);
		END IF;
    END LOOP;
END;
$$;

CREATE EVENT TRIGGER 
	test_event_trigger_for_drops_tr
ON 
	sql_drop
EXECUTE FUNCTION 
	test_event_trigger_for_drops();

CREATE TYPE mmseq_result AS (
	query_id bigint,
	target_id bigint,
	raw_score double precision,
	bit_score double precision,
	e_value double precision,
	q_start integer,
	q_end integer,
	q_len integer,
	t_start integer,
	t_end integer,
	t_len integer,
	q_aln text,
	t_aln text,
	cigar text,
	aln_len integer,
	mismatch integer,
	gap_open integer,
	pident double precision
);

CREATE OR REPLACE FUNCTION nucl_search_one_to_one(
		nucl_seq,
		nucl_seq,
		subst_matrix_name text = 'nucleotide',
		kmer_length integer = 7,
		kmer_gen_threshold integer = 9,
		ungapped_alignment_score integer = 6,
		eval_threshold double precision = 0.001,
		gap_open_cost integer = 5,
		gap_penalty_cost integer = 2,
		thread_number integer = 1
	)
    RETURNS SETOF mmseq_result
    AS 'MODULE_PATHNAME', 'nucl_search_one_to_one' LANGUAGE C;

CREATE OR REPLACE FUNCTION nucl_search_one_to_arr(
		nucl_seq,
		nucl_seq[],
		subst_matrix_name text = 'nucleotide',
		kmer_length integer = 7,
		kmer_gen_threshold integer = 9,
		ungapped_alignment_score integer = 6,
		eval_threshold double precision = 0.001,
		gap_open_cost integer = 5,
		gap_penalty_cost integer = 2,
		thread_number integer = 1
	)
    RETURNS SETOF mmseq_result
    AS 'MODULE_PATHNAME', 'nucl_search_one_to_arr' LANGUAGE C;

CREATE OR REPLACE FUNCTION nucl_search_one_to_db(
		nucl_seq,
		text,
		text,
		target_ids bigint[] = NULL,
		subst_matrix_name text = 'nucleotide',
		kmer_length integer = 7,
		kmer_gen_threshold integer = 9,
		ungapped_alignment_score integer = 6,
		eval_threshold double precision = 0.001,
		gap_open_cost integer = 5,
		gap_penalty_cost integer = 2,
		thread_number integer = 1
	)
    RETURNS SETOF mmseq_result
    AS 'MODULE_PATHNAME', 'nucl_search_one_to_db' LANGUAGE C;

CREATE OR REPLACE FUNCTION nucl_search_arr_to_one(
		nucl_seq[],
		nucl_seq,
		subst_matrix_name text = 'nucleotide',
		kmer_length integer = 7,
		kmer_gen_threshold integer = 9,
		ungapped_alignment_score integer = 6,
		eval_threshold double precision = 0.001,
		gap_open_cost integer = 5,
		gap_penalty_cost integer = 2,
		thread_number integer = 1
	)
    RETURNS SETOF mmseq_result
    AS 'MODULE_PATHNAME', 'nucl_search_arr_to_one' LANGUAGE C;

CREATE OR REPLACE FUNCTION nucl_search_arr_to_arr(
		nucl_seq[],
		nucl_seq[],
		subst_matrix_name text = 'nucleotide',
		kmer_length integer = 7,
		kmer_gen_threshold integer = 9,
		ungapped_alignment_score integer = 6,
		eval_threshold double precision = 0.001,
		gap_open_cost integer = 5,
		gap_penalty_cost integer = 2,
		thread_number integer = 1
	)
    RETURNS SETOF mmseq_result
    AS 'MODULE_PATHNAME', 'nucl_search_arr_to_arr' LANGUAGE C;

CREATE OR REPLACE FUNCTION nucl_search_arr_to_db(
		nucl_seq[],
		text,
		text,
		target_ids bigint[] = NULL,
		subst_matrix_name text = 'nucleotide',
		kmer_length integer = 7,
		kmer_gen_threshold integer = 9,
		ungapped_alignment_score integer = 6,
		eval_threshold double precision = 0.001,
		gap_open_cost integer = 5,
		gap_penalty_cost integer = 2,
		thread_number integer = 1
	)
    RETURNS SETOF mmseq_result
    AS 'MODULE_PATHNAME', 'nucl_search_arr_to_db' LANGUAGE C;

CREATE OR REPLACE FUNCTION nucl_search_db_to_one(
		text,
		text,
		nucl_seq,
		query_ids bigint[] = NULL,
		subst_matrix_name text = 'nucleotide',
		kmer_length integer = 7,
		kmer_gen_threshold integer = 9,
		ungapped_alignment_score integer = 6,
		eval_threshold double precision = 0.001,
		gap_open_cost integer = 5,
		gap_penalty_cost integer = 2,
		thread_number integer = 1
	)
    RETURNS SETOF mmseq_result
    AS 'MODULE_PATHNAME', 'nucl_search_db_to_one' LANGUAGE C;

CREATE OR REPLACE FUNCTION nucl_search_db_to_arr(
		text,
		text,
		nucl_seq[],
		query_ids bigint[] = NULL,
		subst_matrix_name text = 'nucleotide',
		kmer_length integer = 7,
		kmer_gen_threshold integer = 9,
		ungapped_alignment_score integer = 6,
		eval_threshold double precision = 0.001,
		gap_open_cost integer = 5,
		gap_penalty_cost integer = 2,
		thread_number integer = 1
	)
    RETURNS SETOF mmseq_result
    AS 'MODULE_PATHNAME', 'nucl_search_db_to_arr' LANGUAGE C;

CREATE OR REPLACE FUNCTION nucl_search_db_to_db(
		text,
		text,
		text,
		text,
		query_ids bigint[] = NULL,
		target_ids bigint[] = NULL,
		subst_matrix_name text = 'nucleotide',
		kmer_length integer = 7,
		kmer_gen_threshold integer = 9,
		ungapped_alignment_score integer = 6,
		eval_threshold double precision = 0.001,
		gap_open_cost integer = 5,
		gap_penalty_cost integer = 2,
		thread_number integer = 1
	)
    RETURNS SETOF mmseq_result
    AS 'MODULE_PATHNAME', 'nucl_search_db_to_db' LANGUAGE C;

CREATE OR REPLACE FUNCTION aa_search_one_to_one(
		aa_seq,
		aa_seq,
		subst_matrix_name text = 'blosum62',
		kmer_length integer = 7,
		kmer_gen_threshold integer = 19,
		ungapped_alignment_score integer = 15,
		eval_threshold double precision = 0.001,
		gap_open_cost integer = 11,
		gap_penalty_cost integer = 1,
		thread_number integer = 1
	)
    RETURNS SETOF mmseq_result
    AS 'MODULE_PATHNAME', 'aa_search_one_to_one' LANGUAGE C;

CREATE OR REPLACE FUNCTION aa_search_one_to_arr(
		aa_seq,
		aa_seq[],
		subst_matrix_name text = 'blosum62',
		kmer_length integer = 7,
		kmer_gen_threshold integer = 19,
		ungapped_alignment_score integer = 15,
		eval_threshold double precision = 0.001,
		gap_open_cost integer = 11,
		gap_penalty_cost integer = 1,
		thread_number integer = 1
	)
    RETURNS SETOF mmseq_result
    AS 'MODULE_PATHNAME', 'aa_search_one_to_arr' LANGUAGE C;

CREATE OR REPLACE FUNCTION aa_search_one_to_db(
		aa_seq,
		text,
		text,
		target_ids bigint[] = NULL,
		subst_matrix_name text = 'blosum62',
		kmer_length integer = 7,
		kmer_gen_threshold integer = 19,
		ungapped_alignment_score integer = 15,
		eval_threshold double precision = 0.001,
		gap_open_cost integer = 11,
		gap_penalty_cost integer = 1,
		thread_number integer = 1
	)
    RETURNS SETOF mmseq_result
    AS 'MODULE_PATHNAME', 'aa_search_one_to_db' LANGUAGE C;

CREATE OR REPLACE FUNCTION aa_search_arr_to_one(
		aa_seq[],
		aa_seq,
		subst_matrix_name text = 'blosum62',
		kmer_length integer = 7,
		kmer_gen_threshold integer = 19,
		ungapped_alignment_score integer = 15,
		eval_threshold double precision = 0.001,
		gap_open_cost integer = 11,
		gap_penalty_cost integer = 1,
		thread_number integer = 1
	)
    RETURNS SETOF mmseq_result
    AS 'MODULE_PATHNAME', 'aa_search_arr_to_one' LANGUAGE C;

CREATE OR REPLACE FUNCTION aa_search_arr_to_arr(
		aa_seq[],
		aa_seq[],
		subst_matrix_name text = 'blosum62',
		kmer_length integer = 7,
		kmer_gen_threshold integer = 19,
		ungapped_alignment_score integer = 15,
		eval_threshold double precision = 0.001,
		gap_open_cost integer = 11,
		gap_penalty_cost integer = 1,
		thread_number integer = 1
	)
    RETURNS SETOF mmseq_result
    AS 'MODULE_PATHNAME', 'aa_search_arr_to_arr' LANGUAGE C;

CREATE OR REPLACE FUNCTION aa_search_arr_to_db(
		aa_seq[],
		text,
		text,
		target_ids bigint[] = NULL,
		subst_matrix_name text = 'blosum62',
		kmer_length integer = 7,
		kmer_gen_threshold integer = 19,
		ungapped_alignment_score integer = 15,
		eval_threshold double precision = 0.001,
		gap_open_cost integer = 11,
		gap_penalty_cost integer = 1,
		thread_number integer = 1
	)
    RETURNS SETOF mmseq_result
    AS 'MODULE_PATHNAME', 'aa_search_arr_to_db' LANGUAGE C;

CREATE OR REPLACE FUNCTION aa_search_db_to_one(
		text,
		text,
		aa_seq,
		query_ids bigint[] = NULL,
		subst_matrix_name text = 'blosum62',
		kmer_length integer = 7,
		kmer_gen_threshold integer = 19,
		ungapped_alignment_score integer = 15,
		eval_threshold double precision = 0.001,
		gap_open_cost integer = 11,
		gap_penalty_cost integer = 1,
		thread_number integer = 1
	)
    RETURNS SETOF mmseq_result
    AS 'MODULE_PATHNAME', 'aa_search_db_to_one' LANGUAGE C;

CREATE OR REPLACE FUNCTION aa_search_db_to_arr(
		text,
		text,
		aa_seq[],
		query_ids bigint[] = NULL,
		subst_matrix_name text = 'blosum62',
		kmer_length integer = 7,
		kmer_gen_threshold integer = 19,
		ungapped_alignment_score integer = 15,
		eval_threshold double precision = 0.001,
		gap_open_cost integer = 11,
		gap_penalty_cost integer = 1,
		thread_number integer = 1
	)
    RETURNS SETOF mmseq_result
    AS 'MODULE_PATHNAME', 'aa_search_db_to_arr' LANGUAGE C;

CREATE OR REPLACE FUNCTION aa_search_db_to_db(
		text,
		text,
		text,
		text,
		query_ids bigint[] = NULL,
		target_ids bigint[] = NULL,
		subst_matrix_name text = 'blosum62',
		kmer_length integer = 7,
		kmer_gen_threshold integer = 19,
		ungapped_alignment_score integer = 15,
		eval_threshold double precision = 0.001,
		gap_open_cost integer = 11,
		gap_penalty_cost integer = 1,
		thread_number integer = 1
	)
    RETURNS SETOF mmseq_result
    AS 'MODULE_PATHNAME', 'aa_search_db_to_db' LANGUAGE C;

CREATE OR REPLACE FUNCTION amb_aa_search_one_to_one(
		amb_aa_seq,
		amb_aa_seq,
		subst_matrix_name text = 'blosum62',
		kmer_length integer = 7,
		kmer_gen_threshold integer = 19,
		ungapped_alignment_score integer = 15,
		eval_threshold double precision = 0.001,
		gap_open_cost integer = 11,
		gap_penalty_cost integer = 1,
		thread_number integer = 1
	)
    RETURNS SETOF mmseq_result
    AS 'MODULE_PATHNAME', 'amb_aa_search_one_to_one' LANGUAGE C;

CREATE OR REPLACE FUNCTION amb_aa_search_one_to_arr(
		amb_aa_seq,
		amb_aa_seq[],
		subst_matrix_name text = 'blosum62',
		kmer_length integer = 7,
		kmer_gen_threshold integer = 19,
		ungapped_alignment_score integer = 15,
		eval_threshold double precision = 0.001,
		gap_open_cost integer = 11,
		gap_penalty_cost integer = 1,
		thread_number integer = 1
	)
    RETURNS SETOF mmseq_result
    AS 'MODULE_PATHNAME', 'amb_aa_search_one_to_arr' LANGUAGE C;

CREATE OR REPLACE FUNCTION amb_aa_search_one_to_db(
		amb_aa_seq,
		text,
		text,
		target_ids bigint[] = NULL,
		subst_matrix_name text = 'blosum62',
		kmer_length integer = 7,
		kmer_gen_threshold integer = 19,
		ungapped_alignment_score integer = 15,
		eval_threshold double precision = 0.001,
		gap_open_cost integer = 11,
		gap_penalty_cost integer = 1,
		thread_number integer = 1
	)
    RETURNS SETOF mmseq_result
    AS 'MODULE_PATHNAME', 'amb_aa_search_one_to_db' LANGUAGE C;

CREATE OR REPLACE FUNCTION amb_aa_search_arr_to_one(
		amb_aa_seq[],
		amb_aa_seq,
		subst_matrix_name text = 'blosum62',
		kmer_length integer = 7,
		kmer_gen_threshold integer = 19,
		ungapped_alignment_score integer = 15,
		eval_threshold double precision = 0.001,
		gap_open_cost integer = 11,
		gap_penalty_cost integer = 1,
		thread_number integer = 1
	)
    RETURNS SETOF mmseq_result
    AS 'MODULE_PATHNAME', 'amb_aa_search_arr_to_one' LANGUAGE C;

CREATE OR REPLACE FUNCTION amb_aa_search_arr_to_arr(
		amb_aa_seq[],
		amb_aa_seq[],
		subst_matrix_name text = 'blosum62',
		kmer_length integer = 7,
		kmer_gen_threshold integer = 19,
		ungapped_alignment_score integer = 15,
		eval_threshold double precision = 0.001,
		gap_open_cost integer = 11,
		gap_penalty_cost integer = 1,
		thread_number integer = 1
	)
    RETURNS SETOF mmseq_result
    AS 'MODULE_PATHNAME', 'amb_aa_search_arr_to_arr' LANGUAGE C;

CREATE OR REPLACE FUNCTION amb_aa_search_arr_to_db(
		amb_aa_seq[],
		text,
		text,
		target_ids bigint[] = NULL,
		subst_matrix_name text = 'blosum62',
		kmer_length integer = 7,
		kmer_gen_threshold integer = 19,
		ungapped_alignment_score integer = 15,
		eval_threshold double precision = 0.001,
		gap_open_cost integer = 11,
		gap_penalty_cost integer = 1,
		thread_number integer = 1
	)
    RETURNS SETOF mmseq_result
    AS 'MODULE_PATHNAME', 'amb_aa_search_arr_to_db' LANGUAGE C;

CREATE OR REPLACE FUNCTION amb_aa_search_db_to_one(
		text,
		text,
		amb_aa_seq,
		query_ids bigint[] = NULL,
		subst_matrix_name text = 'blosum62',
		kmer_length integer = 7,
		kmer_gen_threshold integer = 19,
		ungapped_alignment_score integer = 15,
		eval_threshold double precision = 0.001,
		gap_open_cost integer = 11,
		gap_penalty_cost integer = 1,
		thread_number integer = 1
	)
    RETURNS SETOF mmseq_result
    AS 'MODULE_PATHNAME', 'amb_aa_search_db_to_one' LANGUAGE C;

CREATE OR REPLACE FUNCTION amb_aa_search_db_to_arr(
		text,
		text,
		amb_aa_seq[],
		query_ids bigint[] = NULL,
		subst_matrix_name text = 'blosum62',
		kmer_length integer = 7,
		kmer_gen_threshold integer = 19,
		ungapped_alignment_score integer = 15,
		eval_threshold double precision = 0.001,
		gap_open_cost integer = 11,
		gap_penalty_cost integer = 1,
		thread_number integer = 1
	)
    RETURNS SETOF mmseq_result
    AS 'MODULE_PATHNAME', 'amb_aa_search_db_to_arr' LANGUAGE C;

CREATE OR REPLACE FUNCTION amb_aa_search_db_to_db(
		text,
		text,
		text,
		text,
		query_ids bigint[] = NULL,
		target_ids bigint[] = NULL,
		subst_matrix_name text = 'blosum62',
		kmer_length integer = 7,
		kmer_gen_threshold integer = 19,
		ungapped_alignment_score integer = 15,
		eval_threshold double precision = 0.001,
		gap_open_cost integer = 11,
		gap_penalty_cost integer = 1,
		thread_number integer = 1
	)
    RETURNS SETOF mmseq_result
    AS 'MODULE_PATHNAME', 'amb_aa_search_db_to_db' LANGUAGE C;
