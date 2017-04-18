require 'rbbt-util'
require 'rbbt/tsv'
require 'rbbt/util/cmd'
require 'rbbt/util/open'
require 'rbbt/resource'
require 'rbbt/sources/organism'

module GDSC
  extend Resource
  self.subdir = "share/databases/GDSC"

  def self.organism
    return Organism.default_code("Hsa")
  end

  GDSC.claim GDSC[".source"].full, :proc do
    CMD.cmd("tr ',' '\\t'|tr '' '\\n'|grep -v '^[[:space:]]'", :in => Open.open("ftp://ftp.sanger.ac.uk/pub4/cancerrxgene/releases/release-5.0/gdsc_manova_input_w5.csv"))
  end

  GDSC.claim GDSC.gene_info, :proc do
    all_fields = TSV.parse_header(GDSC[".source"].full.open, :header_hash => '')
    genes = all_fields.fields[5..71]
    tsv = GDSC[".source"].full.tsv :header_hash => '', :fields => genes, :type => :list
    tsv.to_s.gsub(/^\s+$/,'')
  end

  GDSC.claim GDSC.gene_variants, :proc do
    tsv =  GDSC.gene_info.tsv :fix => proc{|l| l.gsub(/::[^\t]*?(\t|$)/, "\t").gsub(/wt|na/,'')}
    tsv.to_s
  end

  GDSC.claim GDSC.gene_CN, :proc do
    tsv =  GDSC.gene_info.tsv :fix => proc{|l| l.gsub(/\t[^\t]*::/, "\t").gsub(/0<cn<8|nci/, '')}
    ppp tsv

    tsv.to_s
  end

  GDSC.claim GDSC.drug_ic50, :proc do
    all_fields = TSV.parse_header(GDSC[".source"].full.open, :header_hash => '')
    drugs = all_fields.fields[75..212]
    tsv = GDSC[".source"].full.tsv :header_hash => '', :fields => drugs, :type => :list
    tsv.fields = drugs.collect{|d| d.sub('_IC_50','')}
    tsv.transpose("Drug").to_s
  end

  GDSC.claim GDSC.gene_expression, :proc do
    tsv = TSV.open("ftp://ftp.sanger.ac.uk/pub4/cancerrxgene/releases/release-5.0/expU133A.txt.zip", :header_hash => '', :type => :list, :cast => :to_f)
    keys = tsv.keys
    file = Organism.identifiers("Hsa").find.to_s
    counts = TSV.field_match_counts(file, keys)
    probe_id = counts.sort_by{|p| p[1]}.last[0]
    tsv.key_field = probe_id
    tsv.to_s
  end


  #### NEW
  
  require 'rbbt/util/excel2tsv'

  GDSC.claim GDSC.cell_lines, :proc do
    url = "ftp://ftp.sanger.ac.uk/pub/project/cancerrxgene/releases/release-6.0/Cell_Lines_Details.xlsx"
    tsv = TSV.xlsx2tsv(url, :sheet => 0)
    tissue = TSV.xlsx2tsv(url, :sheet => 1)
    tissue.fields = ["COSMIC identifier"] + tissue.fields[2..-1]
    tissue = tissue.reorder "COSMIC identifier"
    tsv = tsv.attach tissue, :persist => false, :persist_data => false
    tsv.to_list.to_s
  end

  GDSC.claim GDSC.drugs, :proc do
    url = "ftp://ftp.sanger.ac.uk/pub/project/cancerrxgene/releases/release-6.0/Screened_Compounds.xlsx"
    tsv = TSV.xlsx2tsv(url, :sheet => 0)
    tsv.to_list.reorder("DRUG NAME").to_s
  end

  GDSC.claim GDSC.drug_equivalences, :proc do
    url = "ftp://ftp.sanger.ac.uk/pub/project/cancerrxgene/releases/release-6.0/GDSC-CCLE-CTRP_conversion.xlsx"
    tsv = TSV.xlsx2tsv(url, :sheet => 1)
    tsv.fields = ["Unknown"] + tsv.fields[2..-1]
    tsv.to_list.to_s
  end

  GDSC.claim GDSC.cell_line_equivalences, :proc do
    url = "ftp://ftp.sanger.ac.uk/pub/project/cancerrxgene/releases/release-6.0/GDSC-CCLE-CTRP_conversion.xlsx"
    tsv = TSV.xlsx2tsv(url, :sheet => 0)
    tsv.to_list.to_s
  end

  GDSC.claim GDSC.gene_CNV, :proc do
    url = "ftp://ftp.sanger.ac.uk/pub/project/cancerrxgene/releases/release-6.0/Gene_level_CN.xlsx"
    tsv = TSV.xlsx2tsv(url, :sheet => 1)
    tsv.delete ""
    tsv.to_list.to_s
  end

  GDSC.claim GDSC.drug_AUC, :proc do
    url = "ftp://ftp.sanger.ac.uk/pub/project/cancerrxgene/releases/release-6.0/v17_fitted_dose_response.xlsx"
    tsv = TSV.xlsx2tsv(url, :sheet => 0, :merge => true, :type => :double)

    tsv.key_field = "ID"
    tsv.fields = %w(IC50_RESULTS_ID COSMIC_ID DRUG_ID MAX_CONC_MICROMOLAR LN_IC50 AUC RMSE)

    tsv.process "DRUG_ID" do |values|
      values.collect do |value|
        value.to_s.sub(/\.0$/,'')
      end
    end

    tsv = tsv.reorder "DRUG_ID", %w(COSMIC_ID MAX_CONC_MICROMOLAR LN_IC50 AUC RMSE), :zipped => true
    tsv.fields = ["COSMIC identifier"] + tsv.fields[2..-1]

    tsv.process "COSMIC identifier" do |value|
      value.to_s.sub(/\.0$/,'')
    end

    tsv.to_s
  end

  GDSC.claim GDSC.cell_line_gene_variants, :proc do
    url = "ftp://ftp.sanger.ac.uk/pub/project/cancerrxgene/releases/release-6.0/WES_variants.xlsx"
    tsv = TSV.xlsx2tsv(url, :sheet => 1, :merge => true)
    tsv.to_s
  end

  GDSC.claim GDSC.cell_line_variants, :proc do

    transcript_cds = Organism.transcript_cds(GDSC.organism).tsv
    transcript_chr = Organism.transcripts(GDSC.organism).tsv :key_field => "Ensembl Transcript ID", :fields => ["Chromosome Name"], :type => :single
    transcript_chr_pos = Organism.transcripts(GDSC.organism).tsv :key_field => "Ensembl Transcript ID", :fields => ["Transcript Start (bp)"], :type => :single, :cast => :to_i
    gene_strand = Organism.gene_positions(GDSC.organism).tsv :key_field => "Ensembl Gene ID", :fields => ["Strand"], :type => :single, :cast => :to_i
    transcript_gene = Organism.transcripts(GDSC.organism).tsv :key_field => "Ensembl Transcript ID", :fields => ["Ensembl Gene ID"], :type => :single
    fields = TSV.parse_header(GDSC.cell_line_gene_variants).fields
    transcript_pos, change_pos = %w(Transcript cDNA).collect{|f| fields.index f}
    dumper = TSV::Dumper.new :key_field => "Sample", :fields => fields +["Genomic Mutations"], :type => :double, :organism => GDSC.organism
    dumper.init
    TSV.traverse GDSC.cell_line_gene_variants, :into => dumper do |sample, values|
      transcripts, changes = values.values_at transcript_pos, change_pos
      mutations = Misc.zip_fields([transcripts, changes]).collect do |transcript, change|
        cds = transcript_cds[transcript]
        chr = transcript_chr[transcript]
        tpos = transcript_chr_pos[transcript]
        gene = transcript_gene[transcript]
        strand = gene_strand[gene]
        mut = Misc.translate_dna_mutation_hgvs2rbbt(change)
        begin
          pos = change.match(/c\.(\d+)/)[1].to_i
          p = strand == 1 ? pos + tpos : tpos - pos
          [chr, p, mut] * ":"
        rescue
          next
        end
      end
      [sample.first, values + [mutations]]
    end

    Misc.collapse_stream(dumper.stream)
  end

end

if __FILE__ == $0
  Log.severity = 0
  Log.tsv GDSC.gene_CN.produce(true).tsv
end
