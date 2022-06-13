CREATE EXTENSION mmseq2;

CREATE TABLE query(
    id bigserial,
    seq amb_aa_seq
);

INSERT INTO query (seq) VALUES ('KDLVKLLKKHGHGPDGPDILTV');

CREATE TABLE target(
    id bigserial,
    seq amb_aa_seq
);

INSERT INTO target (seq) VALUES 
    ('AFSAEDVLKEYDRRRRMEALL'),
    ('SLYYPNDRKLLDYKEWSPPRVQVECPKAPVEWNNP'),
    ('IVGHFSGIKYKGEKAQASEVDVNKMCCWVSKFKDAMR'),
    ('KAYNLTVEGVEGFVRYSRVTKQHVAAFL'),
    ('DFKALVESAHRMRQGHMINVKYILYQ');

SELECT * FROM amb_aa_search_db_to_db('query', 'seq', 'target', 'seq', kmer_gen_threshold => 17);
