'use strict';

var through2 = require('through2');
var ggp = require('gramene-gene-positions');
var _ = require('lodash');

function makeKeys(gene) {
  if (!gene._transcripts) {
    gene._transcripts = _.keyBy(gene.gene_structure.transcripts,'id');
  }
  if (!gene._exons) {
    gene._exons = _.keyBy(gene.gene_structure.exons,'id');
  }
}

function bedGene(gene) {
  return [
    gene.location.region,
    gene.location.start-1,
    gene.location.end,
    gene._id,
    0,
    (gene.location.strand === 1) ? '+' : '-'
  ];
}

function bedTranscript(gene,transcript_id) {
  makeKeys(gene);
  var tr = gene._transcripts[transcript_id];
  var bed = [
    gene.location.region,
    gene.location.start, // might change if alt TSS
    gene.location.end,   // might change if alt TTS
    transcript_id.replace(/__/g,'.'),
    0,
    (gene.location.strand === 1) ? '+' : '-',
    gene.location.start,
    gene.location.end,
    0,
    tr.exons.length,
    '',
    ''
  ];
  var blockSizes = [];
  var blockStarts = [];
  if (gene.location.strand === 1) {
    bed[1] = ggp.remap(gene, 1, 'transcript', 'genome', transcript_id);
    bed[2] = ggp.remap(gene, tr.length, 'transcript', 'genome', transcript_id);
    if (tr.cds) {
      bed[6] = ggp.remap(gene, tr.cds.start, 'transcript', 'genome', transcript_id);
      bed[7] = ggp.remap(gene, tr.cds.end, 'transcript', 'genome', transcript_id);
    }
    else {
      bed[6] = bed[1];
      bed[7] = bed[2];
    }
    // iterate over exons in order
    tr.exons.forEach(function (ex_id) {
      var exon = gene._exons[ex_id];
      blockStarts.push(exon.start-1);
      blockSizes.push(exon.end-exon.start+1);
    });
  }
  else {
    bed[2] = ggp.remap(gene, 1, 'transcript', 'genome', transcript_id);
    bed[1] = ggp.remap(gene, tr.length, 'transcript', 'genome', transcript_id);
    if (tr.cds) {
      bed[7] = ggp.remap(gene, tr.cds.start, 'transcript', 'genome', transcript_id);
      bed[6] = ggp.remap(gene, tr.cds.end, 'transcript', 'genome', transcript_id);
    }
    else {
      bed[6] = bed[1];
      bed[7] = bed[2];
    }
    // iterate over exons in reverse order
    tr.exons.slice().reverse().forEach(function(ex_id) {
      var exon = gene._exons[ex_id];
      blockStarts.push(ggp.remap(gene, exon.end, 'gene', 'genome', transcript_id) - bed[1]);
      blockSizes.push(exon.end-exon.start+1);
    });
  }
  bed[10] = blockSizes.join(',');
  bed[11] = blockStarts.join(',');
  bed[1]--;
  return bed;
}

function bedIntrons(gene, transcript_id) {
  makeKeys(gene);
  var beds=[];
  var tr = gene._transcripts[transcript_id];
  if (gene.location.strand === 1) {
    for(var i=1; i<tr.exons.length; i++) {
      var intron_id = tr.exons[i-1] + '-intron-' + tr.exons[i];
      beds.push([
        gene.location.region,
        ggp.remap(gene, gene._exons[tr.exons[i-1]].end, 'gene', 'genome')-1,
        ggp.remap(gene, gene._exons[tr.exons[i]].start, 'gene', 'genome'),
        intron_id.replace(/__/g,'.'),
        0,
        gene.location.strand
      ]);
    }
  }
  else {
    for(var i=1; i<tr.exons.length; i++) {
      var intron_id = tr.exons[i-1] + '-intron-' + tr.exons[i];
      beds.push([
        gene.location.region,
        ggp.remap(gene, gene._exons[tr.exons[i]].start, 'gene', 'genome')-1,
        ggp.remap(gene, gene._exons[tr.exons[i-1]].end, 'gene', 'genome'),
        intron_id.replace(/__/g,'.'),
        0,
        gene.location.strand
      ]);
    }
  }
  return beds;
}

function bedCDS(gene, transcript_id) {
  makeKeys(gene);
  var beds=[];
  var tr = gene._transcripts[transcript_id];
  if (tr.cds) {
    var cds_start = ggp.remap(gene,tr.cds.start, 'transcript','gene');
    var cds_end = ggp.remap(gene,tr.cds.end, 'transcript','gene');
    tr.exons.forEach(function (exon_id) {
      var exon = gene._exons[exon_id];
      if (exon.end >= cds_start && exon.start <= cds_end) {
        if (gene.location.strand === 1) {
          beds.push([
            gene.location.region,
            (exon.start < cds_start) ?
              ggp.remap(gene, cds_start, 'gene', 'genome')-1
            : ggp.remap(gene, exon.start, 'gene', 'genome')-1,
            (exon.end > cds_end) ?
              ggp.remap(gene, cds_end, 'gene', 'genome')
            : ggp.remap(gene, exon.end, 'gene', 'genome'),
            exon_id.replace(/__/g, '.'),
            0,
            '+'
          ]);
        }
        else {
          beds.push([
            gene.location.region,
            (exon.end > cds_end) ?
              ggp.remap(gene, cds_end, 'gene', 'genome')-1
            : ggp.remap(gene, exon.end, 'gene', 'genome')-1,
            (exon.start < cds_start) ?
              ggp.remap(gene, cds_start, 'gene', 'genome')
            : ggp.remap(gene, exon.start, 'gene', 'genome'),
            exon_id.replace(/__/g, '.'),
            0,
            '-'
          ]);
        }
      }
    });
  }
  return beds;
}

function bedUTR5(gene, transcript_id) {
  makeKeys(gene);
  var beds=[];
  var tr = gene._transcripts[transcript_id];
  if (tr.cds && tr.cds.start > 1) {
    var cds_gene_start = ggp.remap(gene, tr.cds.start, 'transcript', 'gene');
    tr.exons.forEach(function (exon_id) {
      var exon = gene._exons[exon_id];
      if (exon.start < cds_gene_start) {
        if (gene.location.strand === 1) {
          beds.push([
            gene.location.region,
            ggp.remap(gene, exon.start, 'gene', 'genome')-1,
            (exon.end < cds_gene_start) ?
              ggp.remap(gene, exon.end, 'gene', 'genome')
            : ggp.remap(gene, cds_gene_start-1, 'gene', 'genome'),
            exon_id.replace(/__/g, '.'),
            0,
            '+'
          ]);
        }
        else {
          beds.push([
            gene.location.region,
            (exon.end < cds_gene_start) ?
              ggp.remap(gene, exon.end, 'gene', 'genome')-1
            : ggp.remap(gene, cds_gene_start-1, 'gene', 'genome')-1,
            ggp.remap(gene, exon.start, 'gene', 'genome'),
            exon_id.replace(/__/g, '.'),
            0,
            '-'
          ]);
        }
      }
    });
  }
  return beds;
}

function bedUTR3(gene, transcript_id) {
  makeKeys(gene);
  var beds=[];
  var tr = gene._transcripts[transcript_id];
  if (tr.cds && tr.cds.end < tr.length) {
    var cds_gene_end = ggp.remap(gene, tr.cds.end, 'transcript', 'gene');
    tr.exons.forEach(function (exon_id) {
      var exon = gene._exons[exon_id];
      if (exon.end > cds_gene_end) {
        if (gene.location.strand === 1) {
          beds.push([
            gene.location.region,
            (exon.start > cds_gene_end) ?
              ggp.remap(gene, exon.start, 'gene', 'genome')-1
            : ggp.remap(gene, cds_gene_end+1, 'gene', 'genome')-1,
            ggp.remap(gene, exon.end, 'gene', 'genome'),
            exon_id.replace(/__/g, '.'),
            0,
            '+'
          ]);
        }
        else {
          beds.push([
            gene.location.region,
            ggp.remap(gene, exon.end, 'gene', 'genome')-1,
            (exon.start > cds_gene_end) ?
              ggp.remap(gene, exon.start, 'gene', 'genome')
            : ggp.remap(gene, cds_gene_end+1, 'gene', 'genome'),
            exon_id.replace(/__/g, '.'),
            0,
            '-'
          ]);
        }
      }
    });
  }
  return beds;
}

var featureExtractor = {
  gene: function(gene) {
    return [bedGene(gene)];
  },
  transcript: function(gene, mode) {
    var beds = [];
    if (mode === 'canonical') {
      var ct_id = gene.gene_structure.canonical_transcript;
      beds.push(bedTranscript(gene,ct_id));
    }
    else if (mode === 'all') {
      makeKeys(gene);
      for(var t_id in gene._transcripts) {
        beds.push(bedTranscript(gene, t_id));
      }
    }
    else {
      // intersection/union of transcripts not supported
    }
    return beds;
  },
  exon: function(gene, mode) {
    makeKeys(gene);
    var beds = [];
    if (mode === 'canonical') {
      var ct_id = gene.gene_structure.canonical_transcript;
      gene._transcripts[ct_id].exons.forEach(function(exon_id) {
        var exon = gene._exons[exon_id];
        beds.push([
          gene.location.region,
          ggp.remap(gene, exon.start, 'gene', 'genome'),
          ggp.remap(gene, exon.end, 'gene', 'genome'),
          exon_id,
          0,
          gene.location.strand === 1 ? '+' : '-'
        ]);
      });
    }
    else if (mode === 'all') {
      gene.gene_structure.exons.forEach(function(exon) {
        beds.push([
          gene.location.region,
          ggp.remap(gene, exon.start, 'gene', 'genome'),
          ggp.remap(gene, exon.end, 'gene', 'genome'),
          exon.id,
          0,
          gene.location.strand === 1 ? '+' : '-'
        ]);
      });
    }
    else {
      // intersection/union of exons not supported
    }
    return beds;
  },
  intron: function(gene, mode) {
    var beds = [];
    if (mode === 'canonical') {
      var ct_id = gene.gene_structure.canonical_transcript;
      beds = bedIntrons(gene,ct_id);
    }
    else if (mode === 'all') {
      makeKeys(gene);
      for(var t_id in gene._transcripts) {
        bedIntrons(gene, t_id).forEach(function(intron) {
          beds.push(intron);
        });
      }
    }
    else {
      // intersection/union of introns not supported
    }
    return beds;
  },
  CDS: function(gene, mode) {
    var beds = [];
    if (mode === 'canonical') {
      var ct_id = gene.gene_structure.canonical_transcript;
      beds = bedCDS(gene,ct_id);
    }
    else if (mode === 'all') {
      makeKeys(gene);
      for(var t_id in gene._transcripts) {
        bedCDS(gene, t_id).forEach(function(intron) {
          beds.push(intron);
        });
      }
    }
    else {
      // intersection/union of CDS not supported
    }
    return beds;
  },
  utr5: function(gene, mode) {
    var beds = [];
    if (mode === 'canonical') {
      var ct_id = gene.gene_structure.canonical_transcript;
      beds = bedUTR5(gene,ct_id);
    }
    else if (mode === 'all') {
      makeKeys(gene);
      for(var t_id in gene._transcripts) {
        bedUTR5(gene, t_id).forEach(function(utr) {
          beds.push(utr);
        });
      }
    }
    else {
      // intersection/union of utr5 not supported
    }
    return beds;
  },
  utr3: function(gene, mode) {
    var beds = [];
    if (mode === 'canonical') {
      var ct_id = gene.gene_structure.canonical_transcript;
      beds = bedUTR3(gene,ct_id);
    }
    else if (mode === 'all') {
      makeKeys(gene);
      for(var t_id in gene._transcripts) {
        bedUTR3(gene, t_id).forEach(function(utr) {
          beds.push(utr);
        });
      }
    }
    else {
      // intersection/union of utr3 not supported
    }
    return beds;
  }
};

var combiners = {
  canonical: 1,
  all: 1
};

module.exports = function(gene,params) {
  var bedFeature = params.bedFeature || 'gene';
  if (!featureExtractor[bedFeature]) {
    throw new Error('feature "' + bedFeature + '" not valid');
  }
  var extract = featureExtractor[bedFeature];

  var bedCombiner = params.bedCombiner || 'canonical';
  if (!combiners[bedCombiner]) {
    throw new Error('combiner "' + bedCombiner + '" not valid');
  }
  var features = extract(gene, bedCombiner).map(function (feature) {
    return feature.join("\t");
  });
  if (features.length > 0) {
    return features.join("\n") + "\n";
  }
  return '';
}