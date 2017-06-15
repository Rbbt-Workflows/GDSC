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

  require 'rbbt/tsv/excel'

  GDSC.claim GDSC.gene_variants, :proc do
    tsv =  GDSC.gene_info.tsv :fix => proc{|l| l.gsub(/::[^\t]*?(\t|$)/, "\t").gsub(/wt|na/,'')}
    tsv.to_s
  end

  GDSC.claim GDSC.gene_CN, :proc do
    tsv =  GDSC.gene_info.tsv :fix => proc{|l| l.gsub(/\t[^\t]*::/, "\t").gsub(/0<cn<8|nci/, '')}
    tsv.to_s
  end

  GDSC.claim GDSC.cell_lines, :proc do
    url = "ftp://ftp.sanger.ac.uk/pub/project/cancerrxgene/releases/release-6.0/Cell_Lines_Details.xlsx"
    tsv = TSV.xlsx(url, :sheet => 0)
    tsv.fields = ["COSMIC cell line ID"] + tsv.fields[1..-1].collect{|f| f.sub(/\n/,' ') }

    tissue = TSV.xlsx(url, :sheet => 1)
    tissue.fields = ["COSMIC cell line ID"] + tissue.fields[1..-1]
    tissue = tissue.reorder "COSMIC cell line ID"

    tsv = tsv.attach tissue, :persist => false, :persist_data => false
    tsv = tsv.select{|k,v| v.length > 6}
    tsv.to_list.to_s
  end

  GDSC.claim GDSC.cell_line_equivalences, :proc do
    url = "ftp://ftp.sanger.ac.uk/pub/project/cancerrxgene/releases/release-6.0/GDSC-CCLE-CTRP_conversion.xlsx"
    tsv = TSV.xlsx(url, :sheet => 0)
    tsv.key_field = "COSMIC cell line ID"
    tsv.fields = tsv.fields[0..6] + ["CTRP exp id2", "CTRP exp id3"] + [tsv.fields[-1]]
    tsv.to_list.to_s
  end

  GDSC.claim GDSC.drugs, :proc do
    url = "ftp://ftp.sanger.ac.uk/pub/project/cancerrxgene/releases/release-6.0/Screened_Compounds.xlsx"
    tsv = TSV.xlsx(url, :sheet => 0)
    tsv.fields = tsv.fields.collect{|f| f == "TARGET" ? "Target (Associated Gene Name)" : f }
    tsv.to_list.reorder("DRUG NAME").to_s
  end

  # ToDo: Revise fields (extra data)
  GDSC.claim GDSC.drug_equivalences, :proc do
    url = "ftp://ftp.sanger.ac.uk/pub/project/cancerrxgene/releases/release-6.0/GDSC-CCLE-CTRP_conversion.xlsx"
    tsv = TSV.xlsx(url, :sheet => 1)
    tsv.fields = ["Unknown"] + tsv.fields[2..-1]
    tsv.to_list.to_s
  end

  GDSC.claim GDSC.gene_CNV, :proc do
    url = "ftp://ftp.sanger.ac.uk/pub/project/cancerrxgene/releases/release-6.0/Gene_level_CN.xlsx"
    tsv = TSV.xlsx(url, :sheet => 1)
    tsv.delete ""
    tsv.to_list.to_s
  end

  GDSC.claim GDSC.cell_line_gene_variants, :proc do
    url = "ftp://ftp.sanger.ac.uk/pub/project/cancerrxgene/releases/release-6.0/WES_variants.xlsx"
    tsv = TSV.xlsx(url, :sheet => 1, :merge => true)
    tsv.to_s
  end

  GDSC.claim GDSC.cell_line_variants, :proc do

    #transcript_cds = Organism.transcript_cds(GDSC.organism).tsv
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
        #cds = transcript_cds[transcript]
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

  GDSC.claim GDSC.drug_AUC, :proc do
    url = "ftp://ftp.sanger.ac.uk/pub/project/cancerrxgene/releases/release-6.0/v17_fitted_dose_response.xlsx"
    tsv = TSV.xlsx(url, :sheet => 0, :merge => true, :type => :double)

    tsv.key_field = "ID"
    tsv.fields = %w(IC50_RESULTS_ID COSMIC_ID DRUG_ID MAX_CONC_MICROMOLAR LN_IC50 AUC RMSE)

    tsv.process "DRUG_ID" do |values|
      values.collect do |value|
        value.to_s.sub(/\.0$/,'')
      end
    end

    Log.tsv tsv

    tsv = tsv.reorder "DRUG_ID", %w(COSMIC_ID MAX_CONC_MICROMOLAR LN_IC50 AUC RMSE), :zipped => true, :merge => true, :one2one => true
    Log.tsv tsv
    tsv.fields = ["COSMIC cell line ID"] + tsv.fields[1..-1]

    tsv.process "COSMIC cell line ID" do |value|
      value.collect{|v| v.to_s.sub(/\.0$/,'')}
    end

    tsv.to_s
  end

  GDSC.claim GDSC.cell_line_CEFs, :proc do
    url = "http://www.cancerrxgene.org/gdsc1000/GDSC1000_WebResources///Data/BEMs/CellLines/CellLines_Mo_BEMs.zip"
    cell_line_CEFs = TSV.setup({}, :key_field => "COSMIC cell line ID", :fields => ["CEF"], :type => :flat)
    TmpFile.with_file do |directory|
      Path.setup(directory)
      FileUtils.mkdir_p directory unless directory.exists?
      zip_file = directory["file.zip"]
      Open.download(url, zip_file)
      Misc.unzip_in_dir(zip_file, directory.content)
      directory.content.glob("*/*.tsv").each do |file|
        tsv = file.tsv :header_hash => '', :type => :list
        cell_lines = tsv.fields
        tsv.through do |cef, values|
          cell_lines.zip(values).each do |cell_line, status|
            next if status == '0'
            cell_line_CEFs[cell_line] ||= []
            cell_line_CEFs[cell_line] << cef
          end
        end
      end
    end
    cell_line_CEFs.to_s
  end
end

if __FILE__ == $0
  Log.severity = 0
  GDSC.drug_AUC.produce(true) if __FILE__ == $0
  Log.tsv GDSC.gene_CN.produce(true).tsv
end

