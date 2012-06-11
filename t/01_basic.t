use strict;
use warnings;
use Algorithm::FeatureSelection;
use Test::More tests => 7;

my $fs = Algorithm::FeatureSelection->new();
isa_ok( $fs, 'Algorithm::FeatureSelection' );
can_ok($fs, 'new');
can_ok($fs, 'calc_information_gain');
can_ok($fs, 'calc_ig');
can_ok($fs, 'calc_pairwise_mutual_information');
can_ok($fs, 'calc_pmi');
can_ok($fs, 'calc_entropy');

