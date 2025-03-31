-- Create a new table with the correct schema
CREATE TABLE users_mag_new (
    id INTEGER NOT NULL PRIMARY KEY AUTOINCREMENT,
    name varchar(255) NOT NULL,
    taxonomy text NULL,
    sample_name varchar(255) NULL,
    fasta_file text NULL,
    created_at datetime NOT NULL,
    user_id integer NOT NULL REFERENCES auth_user(id) DEFERRABLE INITIALLY DEFERRED,
    fai_file text NULL,
    gff_file text NULL,
    protein_file text NULL,
    quality_data text NULL,
    annotations_data text NULL,
    completeness float NULL,
    contamination float NULL,
    strain_heterogeneity float NULL,
    gene_count integer NULL,
    trna_count integer NULL,
    rrna_count integer NULL,
    cds_count integer NULL,
    annotation_summary text NULL
);

-- Copy data from the old table to the new table
INSERT INTO users_mag_new(
    id, name, taxonomy, sample_name, fasta_file, created_at, user_id, 
    fai_file, gff_file, protein_file, quality_data, annotations_data, 
    completeness, contamination, strain_heterogeneity, gene_count, 
    trna_count, rrna_count, cds_count, annotation_summary
) 
SELECT 
    id, name, taxonomy, sample_name, fasta_file, created_at, user_id, 
    fai_file, gff_file, protein_file, quality_data, annotations_data, 
    completeness, contamination, strain_heterogeneity, gene_count, 
    trna_count, rrna_count, cds_count, annotation_summary
FROM users_mag;

-- Drop the old table
DROP TABLE users_mag;

-- Rename the new table to the original name
ALTER TABLE users_mag_new RENAME TO users_mag;

-- Recreate the index
CREATE INDEX users_mag_user_id_0cdec01e ON users_mag (user_id); 