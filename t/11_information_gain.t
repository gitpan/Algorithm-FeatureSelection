use strict;
use warnings;
use Algorithm::FeatureSelection;
use Test::More tests => 2;

my $fs = Algorithm::FeatureSelection->new();
isa_ok( $fs, 'Algorithm::FeatureSelection' );

my $test_features = {
    'blog_A' => { 'Female' => 3 },
    'blog_B' => { 'Female' => 1 },
    'blog_C' => {
        'Female' => 4,
        'Male'   => 3
    },
    'blog_D' => {
        'Female' => 389,
        'Male'   => 913
    },
    'blog_E' => { 'Male' => 1 },
    'blog_F' => { 'Male' => 1 },
    'blog_G' => {
        'Female' => 3,
        'Male'   => 2
    },
    'blog_H' => { 'Female' => 1 },
    'blog_I' => { 'Female' => 2 },
    'blog_J' => {
        'Female' => 4,
        'Male'   => 10
    },
    'blog_K' => { 'Male' => 1 },
    'blog_L' => {
        'Female' => 1,
        'Male'   => 3
    },
    'blog_M' => { 'Male' => 1 },
    'blog_N' => { 'Male' => 1 },
    'blog_O' => { 'Male' => 1 },
    'blog_P' => {
        'Female' => 7,
        'Male'   => 15
    },
    'blog_Q' => { 'Male'   => 2 },
    'blog_R' => { 'Female' => 1 },
    'blog_S' => { 'Female' => 4 },
    'blog_T' => {
        'Female' => 18,
        'Male'   => 1
    },
    'blog_U' => { 'Male' => 1 },
};

my $correct_ig = {
    'blog_A' => '0.8980',
    'blog_B' => '0.8980',
    'blog_C' => '0.0010',
    'blog_D' => '0.0105',
    'blog_E' => '0.8980',
    'blog_F' => '0.8980',
    'blog_G' => '0.0009',
    'blog_H' => '0.8980',
    'blog_I' => '0.8980',
    'blog_J' => '0.0000',
    'blog_K' => '0.8980',
    'blog_L' => '0.0000',
    'blog_M' => '0.8980',
    'blog_N' => '0.8980',
    'blog_O' => '0.8980',
    'blog_P' => '0.0000',
    'blog_Q' => '0.8980',
    'blog_R' => '0.8980',
    'blog_S' => '0.8980',
    'blog_T' => '0.0182',
    'blog_U' => '0.8980',
};

my $ig = $fs->calc_information_gain($test_features);

for(keys %$ig){
    $ig->{$_} = sprintf("%6.4f", $ig->{$_});
}

is_deeply( $ig, $correct_ig );