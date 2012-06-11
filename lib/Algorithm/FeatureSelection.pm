package Algorithm::FeatureSelection;
use strict;
use warnings;
use List::Util qw(sum);

our $VERSION = '0.01';

sub new {
    my $class = shift;
    my $self = bless {@_}, $class;
    return $self;
}

sub calc_pmi {
    my $self = shift;
    $self->calc_pairewise_mutual_information(@_);
}

sub calc_ig {
    my $self = shift;
    $self->calc_information_gain(@_);
}

sub calc_pairwise_mutual_information {
    my $self     = shift;
    my $features = shift;

    ## -----------------------------------------------------------------
    ##
    ## The argument is expected as below.
    ##
    ## $features = {
    ##     feature_1 => {
    ##         class_a => 10,
    ##         class_b => 2,
    ##     },
    ##     feature_2 => {
    ##         class_b => 11,
    ##         class_d => 32
    ##     },
    ##           .
    ##           .
    ##           .
    ## };
    ##
    ## -----------------------------------------------------------------
    ##
    ## Pairewise Mutual Information
    ##
    ## PMI(w, c) = log ( P( Xw = 1, C = c ) / P( Xw=1 )P( C=c ) )
    ##
    ##    c.f.  w = feature
    ## 		    c = Class
    ##
    ## -----------------------------------------------------------------

    my $feature_count;
    my $class_count;
    my $co_occur_count;
    my $all_features_num;
    while ( my ( $feature, $ref ) = each %$features ) {
        while ( my ( $class, $count ) = each %$ref ) {
            $feature_count->{$feature}                    += $count;
            $class_count->{$class}                        += $count;
            $co_occur_count->{ $class . "\t" . $feature } += $count;
            $all_features_num                             += $count;
        }
    }

    my $PMI;

    for ( keys %$co_occur_count ) {
        my $f12 = $co_occur_count->{$_};
        my ( $class, $feature ) = split "\t", $_;
        my $f1 = $feature_count->{$feature};
        my $f2 = $class_count->{$class};

        my $pmi_score = _log2( ( $f12 / $all_features_num )
            / ( ( $f1 / $all_features_num ) * ( $f2 / $all_features_num ) ) );

        $PMI->{$feature}->{$class} = $pmi_score;
    }

    return $PMI;
}

sub calc_information_gain {
    my $self     = shift;
    my $features = shift;

    ## -----------------------------------------------------------------
    ##
    ## The argument is expected as below.
    ##
    ## $features = {
    ##     feature_1 => {
    ##         class_a => 10,
    ##         class_b => 2,
    ##     },
    ##     feature_2 => {
    ##         class_b => 11,
    ##         class_d => 32
    ##     },
    ##           .
    ##           .
    ##           .
    ## };
    ##
    ## -----------------------------------------------------------------
    ##
    ## Information Gain
    ##
    ## IG(w) = H(C) - ( P(Xw = 1) H(C|Xw = 1) + P(Xw = 0) H(C|Xw = 0) )
    ##
    ##    c.f. ：w = feature
    ## 		     C = class
    ##
    ## -----------------------------------------------------------------

    my $IG;

    my $classes;
    my $classes_sum;
    my $all_features_num;
    while ( my ( $feature, $ref ) = each %$features ) {
        while ( my ( $class, $count ) = each %$ref ) {
            $classes->{$class}->{$feature} += $count;
            $classes_sum->{$class}         += $count;
            $all_features_num              += $count;
        }
    }

    my @array;
    while ( my ( $class, $ref ) = each %$classes ) {
        my $sum     = sum( values %$ref );
        my $p_class = $sum / $all_features_num;
        push @array, $p_class;
    }
    my $entropy = $self->calc_entropy( \@array );

    while ( my ( $feature, $ref ) = each %$features ) {

        my $sum = sum( values %$ref );

        # H ( C | Xw = 1)
        my $on_entropy;
        {
            my @array;
            while ( my ( $class, $count ) = each %$ref ) {
                my $p_class_feature = $count / $sum;
                push @array, $p_class_feature;
            }

            $on_entropy = $self->calc_entropy( \@array );
        }

        # H ( C | Xw = 0)
        my $off_entropy;
        {
            my @array;
            while ( my ( $class, $count ) = each %$ref ) {

                my $p_class_feature = ( $classes_sum->{$class} - $count )
                    / ( $all_features_num - $sum );
                push @array, $p_class_feature;
            }

            $off_entropy = $self->calc_entropy( \@array );
        }

        # Information Gain
        my $ig
            = $entropy
            - ( ( $sum / $all_features_num ) 
            * $on_entropy
                + ( ( $all_features_num - $sum ) / $all_features_num )
                * $off_entropy );

        $IG->{$feature} = $ig;
    }

    return $IG;
}

sub calc_entropy {
    my $self = shift;
    my $data = shift;

    my @ratio;
    if ( ref $data eq 'HASH' ) {
        @ratio = _ratio( [ values %$data ] );
    }
    elsif ( ref $data eq 'ARRAY' ) {
        if ( sum(@$data) == 1 ) {
            @ratio = @$data;
        }
        else {
            @ratio = _ratio($data);
        }
    }

    my $entropy;
    for my $p (@ratio) {
        $entropy += -$p * _log2($p);
    }
    return $entropy;

}

sub _ratio {
    my $arrayref = shift;
    my @ratio;
    my $sum = sum(@$arrayref);
    for (@$arrayref) {
        push @ratio, $_ / $sum;
    }
    return @ratio;
}

sub _log2 {
    my $n = shift;
    log($n) / log(2);
}

1;
__END__

=head1 NAME

Algorithm::FeatureSelection -

=head1 SYNOPSIS

  use Algorithm::FeatureSelection;
  my $fs = Algorithm::FeatureSelection->new();

  # feature-class data structure ...
  my $features = {
    feature_1 => {
        class_a => 10,
        class_b => 2,
    },
    feature_2 => {
        class_b => 11,
        class_d => 32
    },
          .
          .
          .
  };

  # get pairwise-mutula-information
  my $pmi = $fs->calc_pairwise_mutual_information($features);
  my $pmi = $fs->calc_pmi($features); # same above

  # get information-gain 
  my $ig = $fs->calc_information_gain($features);
  my $ig = $fs->calc_ig($features); # same above



=head1 DESCRIPTION

This library is an perl implementation of 'Pairwaise Mutual Information' and 'Information Gain' 
that are used as well-known method of feature selection on text mining fields.

=head1 METHOD

=head2 new()

=head2 calc_information_gain( $features )

  my $features = {
    feature_1 => {
        class_a => 10,
        class_b => 2,
    },
    feature_2 => {
        class_b => 11,
        class_d => 32
    },
          .
          .
          .
  };
  my $fs = Algorithm::FeatureSelection->new();
  my $ig = $fs->information_gain($features);

=head2 calc_ig( $features )

  short name of calc_information_gain()

=head2 calc_pairwise_mutual_information( $features )

  my $features = {
    feature_1 => {
        class_a => 10,
        class_b => 2,
    },
    feature_2 => {
        class_b => 11,
        class_d => 32
    },
          .
          .
          .
  };
  my $fs = Algorithm::FeatureSelection->new();
  my $pmi = $fs->calc_pairwise_mutual_information($features);

=head2 calc_pmi( $features )

  short name of calc_pairwise_mutual_information()

=head2 calc_entropy(HASH|ARRAY)

  calcurate entropy. 

=head1 AUTHOR

Takeshi Miki E<lt>miki@cpan.orgE<gt>

=head1 SEE ALSO

=head1 LICENSE

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself.

=cut
