require 'rbbt-util'
require 'rbbt/tsv'
require 'rbbt/util/cmd'
require 'rbbt/util/open'
require 'rbbt/resource'
require 'rbbt/sources/organism'

module GDSC
  extend Resource
  self.subdir = "share/databases/GDSC"

  GDSC.claim GDSC[".source"].full, :proc do
    CMD.cmd("tr ',' '\\t'|tr '' '\\n'|grep -v '^[[:space:]]'", :in => Open.open("ftp://ftp.sanger.ac.uk/pub4/cancerrxgene/current_release/gdsc_manova_input_w2.csv"))
  end

  GDSC.claim GDSC.gene_info, :proc do
    all_fields = TSV.parse_header(GDSC[".source"].full.open, :header_hash => '')
    genes = all_fields.fields[5..71]
    tsv = GDSC[".source"].full.tsv :header_hash => '', :fields => genes, :type => :list
    tsv.to_s.gsub(/^\s+$/,'')
  end

  GDSC.claim GDSC.gene_variants, :proc do
    tsv =  GDSC.gene_info.tsv :fix => proc{|l| l.gsub(/::.*?(\t|$)/, "\t").gsub(/wt|na/,'')}
    tsv.to_s
  end

  GDSC.claim GDSC.gene_CN, :proc do
    tsv =  GDSC.gene_info.tsv :fix => proc{|l| l.gsub(/\t.*::/, "\t").gsub(/0<cn<8|nci/, '')}
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
    tsv = TSV.open("ftp://ftp.sanger.ac.uk/pub4/cancerrxgene/releases/release-2.0/expU133A.txt", :header_hash => '', :type => :list, :cast => :to_f)
    keys = tsv.keys
    file = Organism.identifiers("Hsa").find.to_s
    counts = TSV.field_match_counts(file, keys)
    probe_id = counts.sort_by{|p| p[1]}.last[0]
    tsv.key_field = probe_id
    tsv.to_s
  end
end

puts GDSC.gene_expression.tsv.summary
