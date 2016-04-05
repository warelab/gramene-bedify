describe('Bedify', function () {
  // json response to http://data.gramene.org/genes?idList=AT1G67500,AT4G25880
  var fixture = require('../support/geneFixture.json');
  var bedify = require('../../index.js');

  var genePlus,geneMinus;

  beforeEach(function() {
    geneMinus = fixture[0];
    genePlus = fixture[1];
  });

  it('should not fail', function () {
    var bedPlus = bedify(genePlus, {bedFeature:'transcript',bedCombiner:'canonical'});
  });
  
  it('should produce 0-based start coords', function () {
    var bedPlus = bedify(genePlus, {bedFeature:'gene'});
    expect(+bedPlus.split("\t")[1]).toEqual(13154707);
  });
  
});