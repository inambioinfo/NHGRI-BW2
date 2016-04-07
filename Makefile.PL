use ExtUtils::MakeMaker;
use strict;
use warnings;

WriteMakefile(
    NAME              => 'NHGRI-BW2',
    VERSION_FROM      => 'lib/NHGRI/BW2/Analysis.pm',
    PREREQ_PM         => {
	},
	BUILD_REQUIRES    => {
	},
    ABSTRACT_FROM     => 'lib/NHGRI/BW2/Analysis.pm',
    AUTHOR            => 'Nancy Fisher Hansen <nhansen@mail.nih.gov>',
    META_MERGE        => {
		resources => {
			repository => {
				type => 'git',
				url  => 'https://github.com/nhansen/NHGRI-BW2.git',
				web  => 'https://github.com/nhansen/NHGRI-BW2',
			},
		},
	},
);

