require 'rbbt/workflow'


require 'rbbt/entity/study'
require 'gdsc'

module GDSC
  extend Workflow

  task :ccle_sample2gdsc => :tsv do
    Workflow.require_workflow "Study"
    s = Study.setup("CCLE")

    gdsc_samples = GDSC.drug_ic52.tsv.fields
    ccle_samples = s.dir.samples.tsv(:fields => []).keys

    ccle_sample2cell_line = s.dir.samples.tsv(:fields => ["Cell line primary name", "Cell line aliases"], :type => :double)

    ccle_sample2gdsc = TSV.setup({}, :key_field => "CCLE Sample", :fields => ["GDSC cell line", "GDSC cell line (soft)"], :type => :single)
    ccle_samples.each do |ccle_sample|
      cell_line_names = ccle_sample2cell_line[ccle_sample].compact.flatten.uniq
      clean_cell_line_names = cell_line_names.collect{|s| s.gsub(/[\s_.-]/,'').downcase  }
      gdsc_sample = gdsc_samples.select{|s| cell_line_names.include? s }
      gdsc_sample_clean = gdsc_samples.select{|s| clean_cell_line_names.include? s.gsub(/[\s_.-]/,'').downcase }
      ccle_sample2gdsc[ccle_sample] = [gdsc_sample.first, gdsc_sample_clean.first]
    end

    ccle_sample2gdsc
  end

  task :gene_samples => :tsv do
    s = Study.setup("CCLE")
    s.gene_sample_matrix.reorder("Ensembl Gene ID")
  end

  dep :ccle_sample2gdsc
  dep :gene_samples
  input :threshold, :float, "Threshold for differential ic50 p-values", 0.1
  task :gene_mutation_sensitivity => :tsv do |threshold|
    ccle_sample2gdsc = step(:ccle_sample2gdsc).load
    gene_samples = step(:gene_samples).load

    samples = gene_samples.fields.sort

    gdsc_samples_in_ccle = ccle_sample2gdsc.values.compact.flatten.uniq

    gene_results = TSV.setup({}, :key_field => "Ensembl Gene ID", :fields => ["Associated Gene Name", "Matching samples", "Growth", "Shrinkage"], :type => :double)

    gene_samples.with_monitor do
      gene_samples.through do |gene, present|
        present_samples = samples.zip(present).select{|s,present| present == "TRUE"}.collect{|sample,p| sample}

        present_gdsc_samples = ccle_sample2gdsc.values_at(*present_samples).compact.flatten.uniq

        rest = gdsc_samples_in_ccle - present_gdsc_samples

        matrix = Matrix.new(GDSC.drug_ic50.produce.find, "Drug")

        diffs = matrix.sample_differences(present_gdsc_samples, rest)

        threshold = threshold.to_f

        growth = Expression.top_up(diffs, threshold).collect{|k,v| [k, v['adjusted.p.values'].to_f]}.sort_by{|p| p.last}.collect{|p| p.first + " (#{p.last})"}
        shrink = Expression.top_down(diffs, threshold).collect{|k,v| [k, v['adjusted.p.values'].to_f]}.sort_by{|p| p.last}.collect{|p| p.first + " (#{p.last})"}

        gene_results[gene] = [gene.name, present_gdsc_samples.length, growth, shrink]
      end
    end

    gene_results
  end

  input :gene, :string, "Gene"
  input :confounding, :array, "Other genes that might be relevant", []
  task :gene_sensitivity => :tsv do |gene, confounding|
    tsv = GDSC.gene_variants.tsv

    all_samples = tsv.keys

    samples = tsv.column(gene).select{|s,v| not v.empty?}.keys

    confounded_samples_tsv = tsv.reorder(:key, confounding)
    confounded_samples = confounded_samples_tsv.select{|s,vs| not vs.flatten.compact.reject{|e| e.empty?}.empty?}.keys

    contrast_samples = all_samples - samples
    case_samples =  samples - confounded_samples

    m = Matrix.new(GDSC.drug_ic50.produce.find, nil, :log, "Drug")
    FileUtils.cp m.differential(case_samples, contrast_samples), self.path
    nil
  end
end

