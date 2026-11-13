#!/usr/bin/perl
# gen.pl: builds the code-derived parts of the nyx C4 set from the Rust source and spec/components.spec.
#
#   perl tools/gen.pl --root <c4 root> --crates <repo>/crates --model <dir for the .dsl fragments> --adr <adr folder relative to the root>
#   perl tools/gen.pl --pages ...    (second phase, after the Structurizr export: builds the pages and the viewer)
#
# Phase 1 (default) writes
#   <model>/<key>.components.dsl   the components of each container (included in the container's block of workspace.dsl)
#   <model>/relationships.dsl      relationships between components, and from components to other containers (read from the source)
#   <model>/views.dsl              one component view per container
#   4-code/<key>/<id>.md           the code page of every component: files, uses, used by, class diagram, public items
#   4-code/<key>/README.md, 4-code/README.md, reports/*.txt
# Phase 2 (--pages) writes 1-context, 2-containers, 3-components, 5-dynamic pages from export/mermaid, and viewer/index.html.
#
# It is a regular-expression reader of Rust (comments, strings and cfg(test) modules blanked, braces counted), not a compiler:
# it reads the public items, their fields, variants and signatures, `impl` blocks, and every `use` and `crate::` / `femcore::` path.
use strict;
use warnings;
use File::Path qw(make_path);
use File::Basename qw(dirname basename);
use Getopt::Long;

my ($root, $crates, $model, $adr_rel, $spec_file, $dsl_file, $pages) = ('.', '../../../crates', 'model', '../adr', 'spec/components.spec', '', 0);
GetOptions('root=s' => \$root, 'crates=s' => \$crates, 'model=s' => \$model, 'adr=s' => \$adr_rel, 'spec=s' => \$spec_file,
           'dsl=s' => \$dsl_file, 'pages' => \$pages) or die "bad options\n";
$spec_file = "$root/$spec_file" unless -f $spec_file;

sub trim { my $s = shift; $s =~ s/^\s+|\s+$//g; return $s }
sub slurp { my $f = shift; open my $h, '<:raw', $f or die "$f: $!"; local $/; my $t = <$h>; close $h; $t =~ s/\r\n/\n/g; return $t }
sub spit { my ($f, $t) = @_; make_path(dirname($f)); open my $h, '>:raw', $f or die "$f: $!"; print $h $t; close $h }

# ---------------------------------------------------------------------------------------------------------- spec
my (@containers, %cont, %crate2key);
{
    open my $fh, '<', $spec_file or die "$spec_file: $!";
    my $cur;
    while (<$fh>) {
        chomp; s/\r$//;
        next if /^\s*#/ || /^\s*$/;
        if (/^\[(.*)\]\s*$/) {
            my ($k, $dir, $crate) = map { trim($_) } split /\|/, $1;
            $cur = { key => $k, dir => $dir, crate => $crate, comps => [], edges => [] };
            push @containers, $cur; $cont{$k} = $cur; $crate2key{$crate} = $k; next;
        }
        if (/^view\s+(.*)$/) { my ($k, $id, $title, $ids) = map { trim($_) } split /\|/, $1; push @{ $cont{$k}{views} }, { id => $id, title => $title, ids => [ split ' ', $ids ] }; next }
        if (/^edge\s+(\S+)\s*->\s*(\S+)\s*\|\s*(.*)$/) { push @{ $cur->{edges} }, { from => $1, to => $2, label => trim($3) }; next }
        my ($id, $title, $files, $adr, $desc) = map { trim($_) } split /\|/, $_, 5;
        my $c = { id => $id, title => $title, files => [ split ' ', $files ], adr => $adr, desc => $desc, key => $cur->{key} };
        push @{ $cur->{comps} }, $c;
    }
}
my %comp;    # "key.id" => comp
for my $k (@containers) { $comp{ "$k->{key}.$_->{id}" } = $_ for @{ $k->{comps} } }

# container descriptions from the DSL (for the pages)
my %cinfo;
if ($dsl_file && -f $dsl_file) {
    my $d = slurp($dsl_file);
    while ($d =~ /^\s*(\w+)\s*=\s*container\s+"([^"]*)"\s+"((?:[^"\\]|\\.)*)"\s+"([^"]*)"/mg) { $cinfo{$1} = { name => $2, desc => $3, tech => $4 } }
}

# ---------------------------------------------------------------------------------------------------------- Rust reading
sub blank_src {
    my $src = shift;
    my $re = qr!(//[^\n]*)|(/\*.*?\*/)|(r(\#*)".*?"\4)|("(?:\\.|[^"\\])*")|('(?:\\(?:u\{[^}]*\}|x[0-9a-fA-F]{2}|.)|[^'\\\n])')!sx;
    (my $c = $src) =~ s{$re}{ my $m = $&; $m =~ s/[^\n]/ /g; $m }gex;
    return $c;
}

sub match_brace {
    my ($r, $i) = @_;
    my $d = 0; pos($$r) = $i;
    while ($$r =~ /([{}])/g) { $d += $1 eq '{' ? 1 : -1; return $-[0] if $d == 0 }
    return length($$r) - 1;
}

# balanced skip of <...> starting at $p (char at $p is '<'); returns the position after the closing '>'
sub skip_angle {
    my ($r, $p) = @_;
    my $d = 0; my $n = length($$r);
    while ($p < $n) {
        my $ch = substr($$r, $p, 1);
        if ($ch eq '-' && substr($$r, $p + 1, 1) eq '>') { $p += 2; next }
        $d++ if $ch eq '<'; $d-- if $ch eq '>';
        $p++;
        return $p if $d == 0;
    }
    return $p;
}

# split on commas outside (), [], {}, <>
sub split_top {
    my $s = shift;
    $s =~ s/->/\x{2192}/g;
    my @out; my $d = 0; my $cur = '';
    for my $ch (split //, $s) {
        if    ($ch =~ /[(\[{<]/) { $d++ }
        elsif ($ch =~ /[)\]}>]/) { $d-- }
        if ($ch eq ',' && $d == 0) { push @out, $cur; $cur = ''; next }
        $cur .= $ch;
    }
    push @out, $cur if $cur =~ /\S/;
    s/\x{2192}/->/g for @out;
    return @out;
}

# a type as a short Mermaid-safe string
sub simple_ty {
    my $t = shift;
    return '' unless defined $t;
    $t =~ s/\s+/ /g; $t = trim($t);
    $t =~ s/^(?:&\s*(?:'\w+\s+)?(?:mut\s+)?)+//;
    $t =~ s/^(?:dyn|impl)\s+//;
    $t =~ s/\s+(?:where|\+).*$//;
    return 'unit' if $t eq '()' || $t eq '';
    return 'tuple' if $t =~ /^\(/;
    return 'slice' if $t =~ /^\[/;
    return 'fn' if $t =~ /^(?:unsafe\s+)?(?:extern\s+\S*\s*)?fn\b/;
    if ($t =~ /^([\w:]+)\s*<\s*(.*)\s*>$/s) {
        my ($b, $in) = ($1, $2); $b =~ s/.*:://;
        my @a = split_top($in);
        @a = grep { !/^'/ } map { trim($_) } @a;
        if (@a == 1 && $a[0] =~ /^[\w:]+$/) { my $x = $a[0]; $x =~ s/.*:://; return "$b~$x~" }
        return $b;
    }
    $t =~ s/.*:://;
    $t =~ s/[^\w~]//g;
    return $t;
}

sub type_words { my $t = shift; return () unless defined $t; return $t =~ /\b([A-Z]\w*)\b/g }

# parse the signature that follows `fn name` at $p in the clean text; returns hash
sub parse_sig {
    my ($r, $p) = @_;
    my $q = $p;
    $q++ while substr($$r, $q, 1) =~ /\s/;
    $q = skip_angle($r, $q) if substr($$r, $q, 1) eq '<';
    $q++ while substr($$r, $q, 1) =~ /\s/;
    my ($params, $ret) = ('', '');
    if (substr($$r, $q, 1) eq '(') {
        my $d = 0; my $s = $q;
        while ($q < length($$r)) {
            my $ch = substr($$r, $q, 1);
            $d++ if $ch eq '('; $d-- if $ch eq ')';
            $q++;
            last if $d == 0;
        }
        $params = substr($$r, $s + 1, $q - $s - 2);
    }
    pos($$r) = $q;
    if ($$r =~ /\G([^{;]*)/g) { my $rest = $1; $rest =~ s/\bwhere\b.*//s; $ret = $1 if $rest =~ /->\s*(.*)$/s }
    my @names;
    for my $pa (split_top($params)) {
        $pa = trim($pa); $pa =~ s/^(?:#\[[^\]]*\]\s*)+//;
        next if $pa =~ /^(?:&\s*(?:'\w+\s+)?)?(?:mut\s+)?self\b/;
        if ($pa =~ /^(?:mut\s+)?(\w+)\s*:/) { push @names, $1 } elsif ($pa =~ /^\(/ || $pa =~ /^&?\w+/) { push @names, '_' }
    }
    return { names => \@names, ret => trim($ret), end => $q };
}

sub parse_file {
    my ($path, $K, $modpath) = @_;
    my $src = slurp($path);
    my @lines = split /\n/, $src, -1;
    my $c = blank_src($src);
    # drop #[cfg(test)] modules
    my $from = 0;
    while ($c =~ /#\[cfg\(test\)\]\s*(?:pub\s+)?mod\s+\w+\s*\{/g) {
        my ($s, $open) = ($-[0], $+[0] - 1);
        my $close = match_brace(\$c, $open);
        my $seg = substr($c, $s, $close - $s + 1); $seg =~ s/[^\n]/ /g; substr($c, $s, $close - $s + 1) = $seg;
        pos($c) = $close + 1;
    }
    my @nl; { my $i = -1; push @nl, $i while ($i = index($c, "\n", $i + 1)) >= 0 }
    my $line_of = sub { my $p = shift; my ($lo, $hi) = (0, scalar @nl); while ($lo < $hi) { my $m = int(($lo + $hi) / 2); if ($nl[$m] < $p) { $lo = $m + 1 } else { $hi = $m } } return $lo };   # 0-based
    my (@bp, @cum); { my $d = 0; while ($c =~ /([{}])/g) { $d += $1 eq '{' ? 1 : -1; push @bp, $-[0]; push @cum, $d } }
    my $depth_at = sub { my $p = shift; my ($lo, $hi) = (0, scalar @bp); while ($lo < $hi) { my $m = int(($lo + $hi) / 2); if ($bp[$m] < $p) { $lo = $m + 1 } else { $hi = $m } } return $lo == 0 ? 0 : $cum[ $lo - 1 ] };
    my $doc_at = sub {
        my $ln = shift; my @d; my $i = $ln - 1;
        $i-- while $i >= 0 && $lines[$i] =~ /^\s*#\[/;
        while ($i >= 0 && $lines[$i] =~ m{^\s*///\s?(.*)$}) { unshift @d, $1; $i-- }
        my ($first) = grep { /\S/ } @d;
        return defined $first ? trim($first) : '';
    };
    my $f = { path => $path, K => $K, mod => $modpath, lines => scalar(@lines), items => [], impls => [], uses => [], clean => $c };
    # header comment
    { my @h; for my $l (@lines) { if ($l =~ m{^\s*//!\s?(.*)$}) { last if trim($1) eq '' && @h; push @h, trim($1) if trim($1) ne ''; last if @h >= 3 } else { last if $l =~ /\S/ || @h } } my $t = join(' ', @h); $t = substr($t, 0, 237) . '...' if length($t) > 240; $f->{header} = $t }
    $f->{header} //= '';

    # impl blocks at depth 0
    my @impl_spans;
    while ($c =~ /^[ \t]*(?:unsafe[ \t]+)?impl\b/mg) {
        my $s = $-[0]; next if $depth_at->($s) != 0;
        my $hs = $+[0]; my $open = index($c, '{', $hs); next if $open < 0;
        my $semi = index($c, ';', $hs); next if $semi >= 0 && $semi < $open;
        my $close = match_brace(\$c, $open);
        my $h = substr($c, $hs, $open - $hs);
        $h =~ s/^\s*<//; if (substr($c, $hs, $open - $hs) =~ /^\s*</) { my $tmp = substr($c, $hs, $open - $hs); my $off = index($tmp, '<'); my $e = skip_angle(\$tmp, $off); $h = substr($tmp, $e) }
        $h =~ s/\bwhere\b.*//s;
        my ($tr, $ty);
        { my $d = 0; my $pos = -1; my $i = 0;
          while ($h =~ /(<|>|\bfor\b)/g) { my $t = $1; if ($t eq '<') { $d++ } elsif ($t eq '>') { $d-- unless substr($h, $-[0] - 1, 1) eq '-' } elsif ($d == 0) { $pos = $-[0]; last } }
          if ($pos >= 0) { $tr = trim(substr($h, 0, $pos)); $ty = trim(substr($h, $pos + 3)) } else { $ty = trim($h) } }
        my $name = simple_ty($ty); $name =~ s/~.*//; $name =~ s/^(?:dyn)//;
        my $trn; if (defined $tr) { $trn = simple_ty($tr); $trn =~ s/~.*//; }
        my $imp = { type => $name, trait => $trn, methods => [], start => $s, end => $close };
        push @impl_spans, $imp;
        push @{ $f->{impls} }, $imp;
    }
    # pub items
    my $item_re = qr/(?:\b(pub(?:\s*\(\s*(?:crate|super|self|in\s+[\w:]+)\s*\))?)\s+|(?<![\w])(?=(?:unsafe\s+|async\s+|const\s+|extern\b[^;{]*?\s)*(?:fn|struct|enum|trait|type|const|static|union)\s))((?:(?:unsafe|async|const|extern(?:\s+\S*)?)\s+)*)(fn|struct|enum|trait|type|const|static|mod|union)\s+(\w+)/;
    while ($c =~ /$item_re/g) {
        my ($vis, $mods, $kind, $name) = (defined $1 ? $1 : '', $2, $3, $4);
        my ($s, $e) = ($-[0], $+[0]);
        next if $kind eq 'mod';
        if ($vis eq '') { my $q = $s - 1; $q-- while $q >= 0 && substr($c, $q, 1) =~ /\s/; my $pc = $q >= 0 ? substr($c, $q, 1) : ';'; next unless $pc =~ /[;{}\]]/ }
        my $dp = $depth_at->($s);
        my $v = $vis eq 'pub' ? '+' : $vis eq '' ? '-' : '~'; next if $name eq '_';
        my $ln = $line_of->($s);
        my $item = { kind => $kind, name => $name, vis => $v, line => $ln + 1, doc => $doc_at->($ln), file => $f };
        if ($dp == 0) {
            if ($kind eq 'fn') { my $sg = parse_sig(\$c, $e); my $raw = substr($c, $s, $sg->{end} - $s); $raw =~ s/\s+/ /g; $item->{sig} = $sg; $item->{raw} = trim($raw . ($sg->{ret} ne '' ? " -> $sg->{ret}" : '')); $item->{raw} = trim(substr($c, $s, 0) . do { my $t = substr($c, $s, $sg->{end} - $s); $t =~ s/\s+/ /g; $t } . ($sg->{ret} ne '' ? " -> $sg->{ret}" : '')) }
            elsif ($kind eq 'struct' || $kind eq 'union') {
                my $q = $e; $q++ while substr($c, $q, 1) =~ /\s/;
                $q = skip_angle(\$c, $q) if substr($c, $q, 1) eq '<';
                $q++ while substr($c, $q, 1) =~ /\s/;
                my $ch = substr($c, $q, 1);
                if ($ch !~ /[{(;]/) { my $o = index($c, '{', $q); my $sm = index($c, ';', $q); $q = ($o >= 0 && ($sm < 0 || $o < $sm)) ? $o : $sm; $ch = substr($c, $q, 1) }
                $item->{fields} = [];
                if ($ch eq '{') {
                    my $cl = match_brace(\$c, $q); my $body = substr($c, $q + 1, $cl - $q - 1);
                    for my $seg (split_top($body)) { $seg =~ s/^(\s*#\[[^\]]*\]\s*)+//; if ($seg =~ /^\s*pub(?:\s*\([^)]*\))?\s+(\w+)\s*:\s*(.+?)\s*$/s) { push @{ $item->{fields} }, { name => $1, type => trim($2) } } }
                } elsif ($ch eq '(') {
                    my $d = 0; my $s2 = $q; while ($q < length($c)) { my $x = substr($c, $q, 1); $d++ if $x eq '('; $d-- if $x eq ')'; $q++; last if $d == 0 }
                    my $i = 0; for my $seg (split_top(substr($c, $s2 + 1, $q - $s2 - 2))) { $seg =~ s/^(\s*#\[[^\]]*\]\s*)+//; if ($seg =~ /^\s*pub(?:\s*\([^)]*\))?\s+(.+?)\s*$/s) { push @{ $item->{fields} }, { name => "_$i", type => trim($1) } } $i++ }
                }
            } elsif ($kind eq 'enum') {
                my $o = index($c, '{', $e); my $cl = match_brace(\$c, $o); my $body = substr($c, $o + 1, $cl - $o - 1);
                $item->{variants} = [];
                for my $seg (split_top($body)) { $seg =~ s/^(\s*#\[[^\]]*\]\s*)+//; if ($seg =~ /^\s*(\w+)\s*(.*)$/s) { push @{ $item->{variants} }, { name => $1, payload => trim($2) } } }
            } elsif ($kind eq 'trait') {
                my $o = index($c, '{', $e); my $cl = match_brace(\$c, $o);
                $item->{methods} = [];
                my $body = substr($c, $o, $cl - $o + 1);
                while ($body =~ /\bfn\s+(\w+)/g) {
                    my $mpos = $o + $-[0]; my $mp = $o + $+[0];
                    next if $depth_at->($mpos) != 1;
                    my $sg = parse_sig(\$c, $mp);
                    push @{ $item->{methods} }, { name => $1, sig => $sg, vis => '+' };
                }
            } elsif ($kind eq 'type' || $kind eq 'const' || $kind eq 'static') {
                my $t = substr($c, $e, 160); $t =~ s/[=;].*//s; $t =~ s/^\s*:\s*//; $item->{ty} = trim($t);
            }
            push @{ $f->{items} }, $item;
        } else {
            # a method of an inherent impl
            for my $imp (@impl_spans) {
                next if defined $imp->{trait};
                if ($s > $imp->{start} && $s < $imp->{end} && $dp == 1 && $kind eq 'fn') {
                    my $sg = parse_sig(\$c, $e); my $raw = substr($c, $s, $sg->{end} - $s); $raw =~ s/\s+/ /g;
                    $item->{sig} = $sg; $item->{raw} = trim($raw . ($sg->{ret} ne '' ? " -> $sg->{ret}" : ''));
                    push @{ $imp->{methods} }, $item; last;
                }
            }
        }
    }
    # use statements (with the enclosing inline module)
    my @inl; while ($c =~ /\bmod\s+(\w+)\s*\{/g) { my $o = $+[0] - 1; push @inl, { name => $1, s => $o, e => match_brace(\$c, $o) } }
    while ($c =~ /(?<![\w:])((?:pub(?:\s*\([^)]*\))?\s+)?)use\s+([^;]+);/g) {
        my ($pubp, $body, $pos) = ($1, $2, $-[0]);
        my $m = $modpath;
        for my $i (sort { $a->{s} <=> $b->{s} } @inl) { $m = ($m eq '' ? '' : "$m\::") . $i->{name} if $pos > $i->{s} && $pos < $i->{e} }
        push @{ $f->{uses} }, { pub => ($pubp =~ /\S/ ? 1 : 0), mod => $m, tree => $body };
    }
    # qualified paths in code
    my %seen;
    while ($c =~ /(?<![\w:])(crate|femcore|geom|nyx_\w+|[a-z_]+)((?:::\w+)+)/g) {
        my ($h, $rest) = ($1, $2);
        next unless $h eq 'crate' || $h eq 'femcore' || $h eq 'geom' || $h =~ /^nyx_/;
        my @s = ($h, grep { length } split /::/, $rest);
        my $k = join '::', @s; next if $seen{$k}++;
        push @{ $f->{paths} }, \@s;
    }
    # a file that exports nothing (a binary's main, a module of private functions) shows its private items; else only the public ones
    unless (grep { $_->{vis} ne '-' } @{ $f->{items} }, map { @{ $_->{methods} } } @{ $f->{impls} }) { }
    else {
        $f->{items} = [ grep { $_->{vis} ne '-' } @{ $f->{items} } ];
        for my $imp (@{ $f->{impls} }) { $imp->{methods} = [ grep { $_->{vis} ne '-' } @{ $imp->{methods} } ] }
    }
    return $f;
}

# ---------------------------------------------------------------------------------------------------------- read all crates
my (%files, %owner, %defs, %reexp, %globs, %modexists);   # files{K}{rel} = parsed file
my %filecomps;     # "K/rel" => [comps]
my @problems;
for my $k (@containers) {
    my $K = $k->{crate}; my $src = "$crates/$k->{dir}/src";
    die "no source folder $src\n" unless -d $src;
    my @rel; { my @st = ('.'); while (@st) { my $d = pop @st; opendir my $dh, "$src/$d" or die; for my $e (sort readdir $dh) { next if $e =~ /^\./; my $p = $d eq '.' ? $e : "$d/$e"; if (-d "$src/$p") { push @st, $p } elsif ($e =~ /\.rs$/) { push @rel, $p } } closedir $dh } }
    @rel = sort @rel;
    for my $rel (@rel) {
        my $mp = $rel; $mp =~ s/\.rs$//; $mp =~ s{/mod$}{}; $mp = '' if $mp eq 'lib'; $mp = 'main' if $mp eq 'main'; $mp =~ s{/}{::}g;
        $mp = '' if $mp eq 'mod';
        # which components
        my @cs;
        for my $c (@{ $k->{comps} }) {
            for my $fs (@{ $c->{files} }) {
                my ($fp) = split /#/, $fs;
                if ($fp =~ m{/$}) { push @cs, $c if index($rel, $fp) == 0 }
                elsif ($fp eq $rel) { push @cs, $c }
            }
        }
        if (!@cs) { push @problems, "unassigned: $K $rel"; next }
        my %u; @cs = grep { !$u{$_->{id}}++ } @cs;
        if (@cs > 1 && !grep { grep { /#/ } @{ $_->{files} } } @cs) { push @problems, "assigned twice: $K $rel to " . join(',', map { $_->{id} } @cs) }
        my $pf = parse_file("$src/$rel", $K, $mp);
        $pf->{rel} = $rel; $pf->{comps} = \@cs; $pf->{sliced} = (@cs > 1) ? 1 : 0;
        $files{$K}{$rel} = $pf; $filecomps{"$K/$rel"} = \@cs;
        $owner{$K}{$mp} //= $cs[0];
        $modexists{$K}{$mp} = 1;
        $defs{$K}{$mp}{ $_->{name} } = 1 for @{ $pf->{items} };
        $defs{$K}{$mp}{ $_->{type} } = 1 for @{ $pf->{impls} };
    }
    # spec entries that match no file
    for my $c (@{ $k->{comps} }) { for my $fs (@{ $c->{files} }) { my ($fp) = split /#/, $fs; my $hit = $fp =~ m{/$} ? scalar(grep { index($_, $fp) == 0 } @rel) : scalar(grep { $_ eq $fp } @rel); push @problems, "spec names no file: $K $fs" unless $hit } }
}
die "component spec does not match the source:\n  " . join("\n  ", @problems) . "\n" if @problems;

# `use` trees -> paths
sub parse_use {
    my ($s) = @_;
    my @tok = $s =~ /(\w+|::|[{},*])/g;
    my $i = 0; my @out; my $rec;
    $rec = sub {
        my @pre = @_; my @segs = @pre;
        while ($i < @tok) {
            my $t = $tok[$i];
            if ($t eq '::') { $i++; next }
            if ($t eq '{') { $i++; while ($i < @tok && $tok[$i] ne '}') { if ($tok[$i] eq ',') { $i++; next } $rec->(@segs) } $i++; return }
            if ($t eq '}' || $t eq ',') { push @out, { path => [@segs] } if @segs > @pre; return }
            if ($t eq '*') { push @out, { path => [ @segs, '*' ] }; $i++; return }
            if ($t eq 'as' && @segs > @pre) { my $al = $tok[ $i + 1 ]; $i += 2; push @out, { path => [@segs], alias => $al }; return }
            push @segs, $t; $i++;
        }
        push @out, { path => [@segs] } if @segs > @pre;
    };
    $rec->();
    for my $o (@out) { pop @{ $o->{path} } if @{ $o->{path} } > 1 && $o->{path}[-1] eq 'self' }
    return @out;
}

# absolute form of a path seen in module $M of crate $K: returns ($K2, [segs]) or ()
sub absolutize {
    my ($K, $M, $segs) = @_;
    my @s = @$segs; return () unless @s;
    my $h = shift @s;
    if ($h eq 'crate') { return ($K, \@s) }
    my @m = $M eq '' ? () : split /::/, $M;
    if ($h eq 'self') { return ($K, [ @m, @s ]) }
    if ($h eq 'super') { pop @m; while (@s && $s[0] eq 'super') { shift @s; pop @m } return ($K, [ @m, @s ]) }
    if (exists $crate2key{$h}) { return ($h, \@s) }
    my $rel = join '::', @m, $h;
    return ($K, [ @m, $h, @s ]) if $modexists{$K}{$rel};
    return ();    # std or an external crate
}

for my $K (keys %files) {
    for my $rel (sort keys %{ $files{$K} }) {
        my $pf = $files{$K}{$rel};
        for my $u (@{ $pf->{uses} }) {
            next unless $u->{pub};
            for my $o (parse_use($u->{tree})) {
                my @p = @{ $o->{path} };
                my $glob = $p[-1] eq '*'; pop @p if $glob;
                next unless @p;
                my ($K2, $abs) = absolutize($K, $u->{mod}, \@p);
                next unless defined $K2;
                if ($glob) { push @{ $globs{$K}{ $u->{mod} } }, join('::', @$abs) if $K2 eq $K; next }
                my $leaf = $o->{alias} // $p[-1];
                $reexp{$K}{ ($u->{mod} eq '' ? '' : "$u->{mod}\::") . $leaf } = { crate => $K2, segs => $abs, home => $pf->{comps}[0] };
            }
        }
    }
}

sub resolve {
    my ($K, $segs, $d) = @_;
    return undef if $d > 10;
    my @s = @$segs;
    for (my $k = @s; $k >= 1; $k--) {
        my $key = join '::', @s[ 0 .. $k - 1 ];
        if (my $r = $reexp{$K}{$key}) {
            return { comp => $r->{home}, mod => '', item => $s[-1], foreign => 1 } if $r->{crate} ne $K;
            my @n = (@{ $r->{segs} }, @s[ $k .. $#s ]);
            return resolve($K, \@n, $d + 1) if join('::', @n) ne join('::', @s);
        }
    }
    for (my $k = @s; $k >= 0; $k--) {
        my $m = join '::', @s[ 0 .. $k - 1 ];
        next unless exists $owner{$K}{$m};
        my $comp = $owner{$K}{$m}; my $item = $k < @s ? $s[$k] : undef;
        if (defined $item && !$defs{$K}{$m}{$item}) {
            for my $g (@{ $globs{$K}{$m} || [] }) { if ($defs{$K}{$g}{$item} && $owner{$K}{$g}) { $comp = $owner{$K}{$g}; $m = $g; last } }
        }
        return { comp => $comp, mod => $m, item => $item };
    }
    return undef;
}

# ---------------------------------------------------------------------------------------------------------- relationships
my %rel;        # "src.id" => { "dst" => { mods=>{}, items=>{} } } ; dst is "key.id" (same container) or "@key" (other container)
my %used_by;    # "key.id" => { "src key.id" => {mods,items} }
my %crossk;     # container-level edges found: "a->b"
sub note_ref {
    my ($srcc, $K, $segs) = @_;
    my $r = resolve($K, $segs, 0);
    my $tk = $crate2key{$K};
    my ($dstkey, $dstcomp, $mod, $item);
    if ($r) { $dstcomp = $r->{comp}; $mod = $r->{mod}; $item = $r->{item} }
    $item = undef if defined $item && $item =~ /^[a-z_0-9]+$/ && $modexists{$K}{ join('::', grep { length } ($mod // ''), $item) };
    my $skey = $srcc->{key};
    if ($tk ne $skey) {
        $crossk{"$skey->$tk"}++;
        my $dst = '@' . $tk;
        my $m = $r ? ($mod eq '' ? 'crate root' : $mod) : 'crate root';
        $rel{"$skey.$srcc->{id}"}{$dst}{mods}{$m}++;
        $rel{"$skey.$srcc->{id}"}{$dst}{items}{$item}++ if defined $item;
        if ($dstcomp) { $used_by{"$tk.$dstcomp->{id}"}{"$skey.$srcc->{id}"}{mods}{$m}++; $used_by{"$tk.$dstcomp->{id}"}{"$skey.$srcc->{id}"}{items}{$item}++ if defined $item }
        else { $used_by{"\@$tk"}{"$skey.$srcc->{id}"}{mods}{$m}++ }
    } elsif ($dstcomp && $dstcomp->{id} ne $srcc->{id}) {
        my $dst = "$tk.$dstcomp->{id}";
        $rel{"$skey.$srcc->{id}"}{$dst}{mods}{$mod eq '' ? 'crate root' : $mod}++;
        $rel{"$skey.$srcc->{id}"}{$dst}{items}{$item}++ if defined $item;
        $used_by{$dst}{"$skey.$srcc->{id}"}{mods}{$mod eq '' ? 'crate root' : $mod}++;
        $used_by{$dst}{"$skey.$srcc->{id}"}{items}{$item}++ if defined $item;
    }
}
for my $K (sort keys %files) {
    for my $rel (sort keys %{ $files{$K} }) {
        my $pf = $files{$K}{$rel};
        next if $pf->{sliced};
        my $srcc = $pf->{comps}[0];
        for my $u (@{ $pf->{uses} }) {
            for my $o (parse_use($u->{tree})) {
                my @p = @{ $o->{path} }; pop @p if $p[-1] eq '*';
                next unless @p;
                my ($K2, $abs) = absolutize($K, $u->{mod}, \@p);
                next unless defined $K2;
                note_ref($srcc, $K2, $abs);
            }
        }
        for my $p (@{ $pf->{paths} || [] }) {
            my ($K2, $abs) = absolutize($K, $pf->{mod}, $p);
            next unless defined $K2;
            note_ref($srcc, $K2, $abs);
        }
    }
}
# hand-written edges of the spec
my %manual;
for my $k (@containers) {
    for my $e (@{ $k->{edges} }) {
        my $from = $e->{from} =~ /\./ ? $e->{from} : "$k->{key}.$e->{from}";
        $manual{$from}{ $e->{to} } = $e->{label};
    }
}

# ---------------------------------------------------------------------------------------------------------- phase 2: pages
if ($pages) { phase_pages(); exit 0 }

# ---------------------------------------------------------------------------------------------------------- DSL output
sub dq { my $s = shift; $s =~ s/\\/\\\\/g; $s =~ s/"/\\"/g; return $s }
sub listshort { my ($h, $n) = @_; my @k = sort keys %$h; my $more = @k > $n ? ' +' . (@k - $n) : ''; @k = @k[ 0 .. $n - 1 ] if @k > $n; return join(', ', @k) . $more }
my $total_loc = 0; my $ncomp = 0;
for my $k (@containers) {
    my $out = "// generated by tools/gen.pl from spec/components.spec and the Rust source; do not edit\n";
    for my $c (@{ $k->{comps} }) {
        my @fl = file_list($k, $c);
        my $loc = 0; $loc += $_->{lines} for @fl;
        $c->{loc} = $loc; $total_loc += $loc; $ncomp++;
        my $code = join(', ', map { my $r = $_->{rel}; "crates/$k->{dir}/src/$r" } @fl);
        $code = shorten_code($k, $c, \@fl);
        $out .= "$c->{id} = component \"" . dq($c->{title}) . "\" \"" . dq($c->{desc}) . "\" \"Rust modules\" {\n";
        $out .= "    properties {\n        \"adr\" \"" . dq($c->{adr}) . "\"\n        \"code\" \"" . dq($code) . "\"\n        \"lines\" \"$loc\"\n        \"codeDiagram\" \"4-code/$k->{key}/$c->{id}.md\"\n    }\n}\n";
    }
    spit("$model/$k->{key}.components.dsl", $out);
}
sub file_list { my ($k, $c) = @_; my @o; for my $rel (sort keys %{ $files{ $k->{crate} } }) { my $pf = $files{ $k->{crate} }{$rel}; push @o, $pf if grep { $_ == $c } @{ $pf->{comps} } } return @o }
sub shorten_code {
    my ($k, $c, $fl) = @_;
    my @r = map { $_->{rel} } @$fl;
    my $s = "crates/$k->{dir}/src: " . join(', ', @r);
    $s = "crates/$k->{dir}/src: " . join(', ', @r[ 0 .. 5 ]) . ", ... (" . scalar(@r) . " files)" if @r > 8;
    return $s;
}

my $rels = "// generated by tools/gen.pl from the use and path references of the Rust source and the edges of spec/components.spec\n";
my %emitted;
my @cross_comp_edges;
for my $src (sort keys %rel) {
    my ($sk) = split /\./, $src;
    for my $dst (sort keys %{ $rel{$src} }) {
        my $h = $rel{$src}{$dst};
        my $mods = listshort($h->{mods}, 2);
        my $items = '';
        my ($target, $label, $tech);
        if ($dst =~ /^\@(\w+)$/) { $target = "nyx.$1"; $label = "uses $1: $mods$items"; $tech = "Rust paths"; }
        else { $target = "nyx.$dst"; $label = "uses $mods$items"; $tech = "Rust paths" }
        $label = substr($label, 0, 57) . '...' if length($label) > 60;
        $rels .= "nyx.$src -> $target \"" . dq($label) . "\" \"$tech\"\n";
        $emitted{"$src->$dst"} = 1;
    }
}
for my $src (sort keys %manual) {
    for my $dst (sort keys %{ $manual{$src} }) {
        my $d = $dst; my $tc; my ($sk) = split /\./, $src;
        if ($d =~ /^\@(\w+)$/) { $tc = "nyx.$1" }
        else { my ($dk) = split /\./, $d; $d = "\@$dk" if $dk ne $sk && $d =~ /\./ ; $tc = $d =~ /^\@(\w+)$/ ? "nyx.$1" : "nyx.$d" }
        next if $emitted{"$src->$d"};
        $rels .= "nyx.$src -> $tc \"" . dq($manual{$src}{$dst}) . "\" \"declared in spec\"\n";
        $rel{$src}{$d}{declared} = $manual{$src}{$dst};
    }
}
spit("$model/relationships.dsl", $rels);

my $views = "// generated by tools/gen.pl\n";
for my $k (@containers) {
    my $name = $cinfo{ $k->{key} } ? $cinfo{ $k->{key} }{name} : $k->{key};
    $views .= "component nyx.$k->{key} \"Components_$k->{key}\" \"The components of $name and what they use.\" {\n    include *\n    autolayout tb\n}\n";
    for my $v (@{ $k->{views} || [] }) {
        $views .= "component nyx.$k->{key} \"Components_$k->{key}_$v->{id}\" \"$v->{title}: a part of the components of $name.\" {\n" . join('', map { "    include nyx.$k->{key}.$_\n" } @{ $v->{ids} }) . "    autolayout tb\n}\n";
    }
}
spit("$model/views.dsl", $views);

# ---------------------------------------------------------------------------------------------------------- code pages
sub adr_link {
    my $n = shift; my $f = sprintf('%04d', $n);
    my ($m) = glob("$root/$adr_rel/$f-*.md");
    return "ADR $n" unless $m;
    return "[ADR $n](../../$adr_rel/" . basename($m) . ")";
}
sub mm_id { my $s = shift; $s =~ s/[^A-Za-z0-9_]/_/g; return $s }
sub md_esc {
    my $s = shift; my @p = split /(`[^`]*`)/, $s; my $o = '';
    for my $x (@p) { if ($x =~ /^`/) { $x =~ s/\|/\\|/g } else { $x =~ s/\|/\\|/g; $x =~ s/\]\(/] (/g; $x =~ s/</&lt;/g; $x =~ s/>/&gt;/g } $o .= $x }
    return $o;
}

my %code_diagrams;    # "key.id" => mermaid text (for the viewer)
my @code_index;
my ($n_classes, $n_items) = (0, 0);
for my $k (@containers) {
    my $K = $k->{crate};
    for my $c (@{ $k->{comps} }) {
        my @fl = file_list($k, $c);
        my %slice;   # for sliced files
        for my $fs (@{ $c->{files} }) { my ($fp, $names) = split /#/, $fs; if (defined $names) { $slice{$fp} = { map { $_ => 1 } split /,/, $names } } }
        # collect
        my (%types, @type_order, %mod_fns, @pub_rows);
        my $ntypes_total = 0;
        for my $pf (@fl) {
            my $sl = $slice{ $pf->{rel} };
            for my $it (@{ $pf->{items} }) {
                next if $sl && !$sl->{ $it->{name} };
                if ($it->{kind} =~ /^(struct|enum|trait|union)$/) {
                    my $t = $types{ $it->{name} } ||= { name => $it->{name}, kind => $it->{kind}, item => $it, methods => [], file => $pf->{rel}, impls => [] };
                    push @type_order, $t->{name} if @{ $t->{methods} } == 0 && !$t->{ordered}++;
                    $t->{item} = $it;
                } else { push @{ $mod_fns{ $pf->{rel} } }, $it }
                push @pub_rows, { it => $it, pf => $pf };
            }
            for my $imp (@{ $pf->{impls} }) {
                my $t = $types{ $imp->{type} };
                if (!$t) { $t = $types{ $imp->{type} } = { name => $imp->{type}, kind => 'struct', item => undef, methods => [], file => $pf->{rel}, impls => [], stub => 1 }; push @type_order, $imp->{type} if !$t->{ordered}++ }
                if (defined $imp->{trait}) { push @{ $t->{impls} }, $imp->{trait}; next }
                for my $m (@{ $imp->{methods} }) {
                    next if $sl && !$sl->{"$imp->{type}::$m->{name}"};
                    push @{ $t->{methods} }, $m; push @pub_rows, { it => $m, pf => $pf, owner => $imp->{type} };
                }
            }
        }
        # drop stubs that have no methods and no traits in this component
        for my $n (keys %types) { delete $types{$n} if $types{$n}{stub} && !@{ $types{$n}{methods} } && !@{ $types{$n}{impls} } }
        # trait methods go to their traits
        # ---- diagram
        my @classes;
        my $multi = @fl > 1;
        my $cap = 36;
        my @ordered = grep { $types{$_} } @type_order; my %seen; @ordered = grep { !$seen{$_}++ } @ordered;
        # classes sorted by file then name, keep the first $cap
        my @cl_names = @ordered;
        my $shown_types = @cl_names > $cap ? $cap : scalar @cl_names;
        my $trunc = @cl_names > $cap ? @cl_names - $cap : 0;
        @cl_names = @cl_names[ 0 .. $shown_types - 1 ];
        my %inshow = map { $_ => 1 } @cl_names;
        my $d = "classDiagram\n";
        my %byfile;
        for my $n (@cl_names) { push @{ $byfile{ $types{$n}{file} } }, $n }
        for my $rel (sort keys %mod_fns) { $byfile{$rel} ||= [] }
        my @rels_lines;
        my $emit_class = sub {
            my ($n, $ind) = @_; my $t = $types{$n}; my $it = $t->{item};
            my $o = "$ind" . "class " . mm_id($n) . " {\n";
            my $kindname = $t->{stub} ? 'impl block' : $t->{kind} eq 'trait' ? 'interface' : $t->{kind} eq 'enum' ? 'enumeration' : 'struct';
            $o .= "$ind    <<$kindname>>\n";
            my $cnt = 0; my $max = 8;
            if ($it && $it->{kind} eq 'enum') {
                my @v = @{ $it->{variants} || [] }; my $shown = 0;
                for my $v (@v) { last if $shown >= 12; $o .= "$ind    $v->{name}\n"; $shown++ }
                $o .= "$ind    more_" . (@v - $shown) . "_variants\n" if @v > $shown;
                for my $v (@v) { for my $w (type_words($v->{payload})) { push @rels_lines, [ $n, $w, 'variant' ] } }
            } elsif ($it && $it->{kind} eq 'trait') {
                my @m = @{ $it->{methods} || [] }; my $shown = 0;
                for my $m (@m) { last if $shown >= $max; $o .= "$ind    +" . mm_id($m->{name}) . "(" . join(', ', @{ $m->{sig}{names} }) . ")" . ($m->{sig}{ret} ne '' ? ' ' . simple_ty($m->{sig}{ret}) : '') . "\n"; $shown++ }
                $o .= "$ind    +more_" . (@m - $shown) . "_methods()\n" if @m > $shown;
            } elsif ($it) {
                my @f = @{ $it->{fields} || [] }; my $shown = 0;
                for my $f (@f) { last if $shown >= $max; $o .= "$ind    +" . mm_id($f->{name}) . ' ' . simple_ty($f->{type}) . "\n"; $shown++; for my $w (type_words($f->{type})) { push @rels_lines, [ $n, $w, 'field' ] } }
                $o .= "$ind    +more_" . (@f - $shown) . "_fields\n" if @f > $shown;
            }
            if ($t->{kind} ne 'trait') {
                my @m = @{ $t->{methods} }; my $shown = 0;
                for my $m (@m) { last if $shown >= $max; my $sg = $m->{sig}; $o .= "$ind    " . ($m->{vis}) . mm_id($m->{name}) . "(" . join(', ', @{ $sg->{names} }) . ")" . ($sg->{ret} ne '' ? ' ' . simple_ty($sg->{ret}) : '') . "\n"; $shown++ }
                $o .= "$ind    +more_" . (@m - $shown) . "_methods()\n" if @m > $shown;
                for my $tr (@{ $t->{impls} }) { push @rels_lines, [ $n, $tr, 'impl' ] }
            }
            $o .= "$ind}\n";
            return $o;
        };
        for my $rel (sort keys %byfile) {
            my $ind = $multi ? '    ' : '';
            my $ns = mm_id($rel);
            $d .= "  namespace $ns {\n" if $multi;
            my $ii = $multi ? '    ' : '  ';
            for my $n (@{ $byfile{$rel} }) { $d .= $emit_class->($n, $ii) }
            if ($mod_fns{$rel} && @{ $mod_fns{$rel} }) {
                my @fns = grep { $_->{kind} eq 'fn' } @{ $mod_fns{$rel} };
                my @oth = grep { $_->{kind} ne 'fn' } @{ $mod_fns{$rel} };
                my $o = "$ii" . "class fns_" . $ns . " {\n$ii    <<module " . $rel . ">>\n";
                my $shown = 0;
                for my $f (@fns) { last if $shown >= 10; $o .= "$ii    " . $f->{vis} . mm_id($f->{name}) . "(" . join(', ', @{ $f->{sig}{names} }) . ")" . ($f->{sig}{ret} ne '' ? ' ' . simple_ty($f->{sig}{ret}) : '') . "\n"; $shown++ }
                $o .= "$ii    +more_" . (@fns - $shown) . "_functions()\n" if @fns > $shown;
                my $sh2 = 0;
                for my $f (@oth) { last if $sh2 >= 4; $o .= "$ii    " . $f->{vis} . mm_id($f->{name}) . ' ' . $f->{kind} . "\n"; $sh2++ }
                $o .= "$ii    +more_" . (@oth - $sh2) . "_items\n" if @oth > $sh2;
                $o .= "$ii}\n";
                $d .= $o if @fns || @oth;
            }
            $d .= "  }\n" if $multi;
        }
        my %relseen;
        for my $r (@rels_lines) {
            my ($a, $b, $kind) = @$r;
            next if $a eq $b; next unless $inshow{$b} || ($kind eq 'impl' && $inshow{$b});
            next unless $inshow{$a};
            my $key = "$a>$b>$kind"; next if $relseen{$key}++;
            if ($kind eq 'impl') { $d .= "  " . mm_id($a) . " ..|> " . mm_id($b) . "\n" }
            else { $d .= "  " . mm_id($a) . " --> " . mm_id($b) . "\n" }
        }
        # implemented traits that are not classes in this diagram: a note per class would be noisy; they are in the table
        my $has_body = ($d =~ tr/\n//) > 1 && $d =~ /class /;
        $d = "classDiagram\n  class no_public_items\n" unless $has_body;
        $n_classes += scalar(@cl_names);
        $code_diagrams{"$k->{key}.$c->{id}"} = $d;

        # ---- page
        my @adr = map { adr_link($_) } grep { /^\d+$/ } map { trim($_) } split /,/, $c->{adr};
        my $p = "# $c->{title}: code\n\n";
        my $cname = $cinfo{ $k->{key} } ? $cinfo{ $k->{key} }{name} : $K;
        $p .= "C4 level 4 (code) of the component **$c->{title}** (`$k->{key}.$c->{id}`) in the container **$cname** (crate `$K`). [Component view of the container](../../3-components/$k->{key}.md) · [All code pages](../README.md)\n\n";
        $p .= "> $c->{desc}\n\n";
        $p .= "**Decisions:** " . (@adr ? join(', ', @adr) : 'none recorded') . "  \n";
        $p .= "**Size:** " . scalar(@fl) . " file" . (@fl == 1 ? '' : 's') . ", $c->{loc} lines.\n\n";
        $p .= "## Files\n\n| File | Lines | What it says about itself |\n|---|---:|---|\n";
        for my $pf (@fl) { $p .= "| `crates/$k->{dir}/src/$pf->{rel}` | $pf->{lines} | " . md_esc($pf->{header}) . " |\n" }
        # uses
        my @uses = sort keys %{ $rel{"$k->{key}.$c->{id}"} || {} };
        $p .= "\n## What it uses\n\n";
        if (@uses) {
            $p .= "| Target | Modules | Items |\n|---|---|---|\n";
            for my $u (@uses) {
                my $h = $rel{"$k->{key}.$c->{id}"}{$u};
                my ($tn, $tl);
                if ($u =~ /^\@(\w+)$/) { my $dk = $1; $tn = "container `$dk`" . ($cinfo{$dk} ? " ($cinfo{$dk}{name})" : ''); $tl = (-f "$root/3-components/$dk.md" || $cont{$dk}) ? "[$tn](../../3-components/$dk.md)" : $tn }
                else { my ($dk, $di) = split /\./, $u, 2; my $dc = $comp{$u}; $tl = $dc ? ($dk eq $k->{key} ? "[$dc->{title}](./$di.md)" : "[$dk: $dc->{title}](../$dk/$di.md)") : $u }
                $p .= "| $tl | " . md_esc(listshort($h->{mods} || {}, 8)) . ($h->{declared} ? " (declared: $h->{declared})" : '') . " | " . md_esc(listshort($h->{items} || {}, 10)) . " |\n";
            }
        } else { $p .= "Nothing outside the component (the standard library and external crates are not listed).\n" }
        my @ub = sort keys %{ $used_by{"$k->{key}.$c->{id}"} || {} };
        $p .= "\n## Used by\n\n";
        if (@ub) {
            $p .= "| Component | Modules | Items |\n|---|---|---|\n";
            for my $u (@ub) {
                my $h = $used_by{"$k->{key}.$c->{id}"}{$u}; my ($sk, $si) = split /\./, $u, 2; my $sc = $comp{$u};
                my $tl = $sk eq $k->{key} ? "[$sc->{title}](./$si.md)" : "[$sk: $sc->{title}](../$sk/$si.md)";
                $p .= "| $tl | " . md_esc(listshort($h->{mods}, 6)) . " | " . md_esc(listshort($h->{items} || {}, 8)) . " |\n";
            }
        } else { $p .= "Nothing else in the workspace names this component (entry points are reached by the program or by tests).\n" }
        $p .= "\n## Class diagram\n\nPublic types with their public fields (or variants) and methods, free functions per file as a `module` box; arrows are field types and `impl Trait`. Capped at $cap types and a few members each; the table below is complete.";
        $p .= " **This component has $trunc more types than the diagram shows.**" if $trunc;
        $p .= "\n\n```mermaid\n$d```\n";
        # items
        $p .= "\n## Public items\n\n| Kind | Item | File:line | Signature or summary |\n|---|---|---|---|\n";
        my $rows = 0; my $cap_rows = 400;
        for my $r (@pub_rows) {
            my $it = $r->{it}; $n_items++;
            last if ++$rows > $cap_rows;
            my $name = defined $r->{owner} ? "$r->{owner}::$it->{name}" : $it->{name};
            my $sum;
            if ($it->{kind} eq 'fn') { $sum = "`" . $it->{raw} . "`" }
            elsif ($it->{kind} eq 'struct') { my @f = @{ $it->{fields} || [] }; $sum = @f ? 'fields: ' . join(', ', map { "`$_->{name}`" } @f[ 0 .. ($#f > 11 ? 11 : $#f) ]) . (@f > 12 ? ', ...' : '') : '' }
            elsif ($it->{kind} eq 'enum') { my @v = map { $_->{name} } @{ $it->{variants} || [] }; $sum = @v ? 'variants: ' . join(', ', map { "`$_`" } @v[ 0 .. ($#v > 11 ? 11 : $#v) ]) . (@v > 12 ? ', ...' : '') : '' }
            elsif ($it->{kind} eq 'trait') { my @m = map { $_->{name} } @{ $it->{methods} || [] }; $sum = @m ? 'methods: ' . join(', ', map { "`$_`" } @m[ 0 .. ($#m > 11 ? 11 : $#m) ]) . (@m > 12 ? ', ...' : '') : '' }
            else { $sum = $it->{ty} ? "`" . $it->{ty} . "`" : '' }
            my $doc = $it->{doc} ne '' ? ' ' . $it->{doc} : '';
            $sum .= ($sum ne '' && $doc ne '' ? ' — ' : '') . $doc;
            $p .= "| $it->{kind} | `" . md_esc($name) . "` | `" . $r->{pf}{rel} . ":$it->{line}` | " . md_esc($sum) . " |\n";
        }
        $p .= "\n(" . (scalar(@pub_rows) - $cap_rows) . " more items not listed.)\n" if @pub_rows > $cap_rows;
        $p .= "\nTrait implementations in this component: " . join(', ', map { "`$_->{type}: $_->{trait}`" } grep { defined $_->{trait} } map { @{ $_->{impls} } } @fl) . "\n" if grep { defined $_->{trait} } map { @{ $_->{impls} } } @fl;
        spit("$root/4-code/$k->{key}/$c->{id}.md", $p);
        push @code_index, [ $k, $c, scalar(@fl), scalar(@pub_rows), scalar(@cl_names), $trunc ];
    }
}
# indexes
for my $k (@containers) {
    my $name = $cinfo{ $k->{key} } ? $cinfo{ $k->{key} }{name} : $k->{crate};
    my $p = "# Code pages of $name\n\nOne page per component, C4 level 4. [Component view](../../3-components/$k->{key}.md) · [Index](../README.md)\n\n| Component | Files | Lines | Public items | Types in the diagram |\n|---|---:|---:|---:|---:|\n";
    for my $r (grep { $_->[0] == $k } @code_index) { my ($kk, $c, $nf, $ni, $nc, $tr) = @$r; $p .= "| [$c->{title}](./$c->{id}.md) | $nf | $c->{loc} | $ni | $nc" . ($tr ? " (+$tr not shown)" : '') . " |\n" }
    spit("$root/4-code/$k->{key}/README.md", $p);
}
{
    my $p = "# Code level (C4 level 4)\n\nOne page per component of every container: the files, what the component uses and what uses it (read from the `use` and path references of the source), a class diagram of its public types and a complete table of its public items. Generated by `tools/gen.pl` from the source; do not edit.\n\n";
    for my $k (@containers) { my @r = grep { $_->[0] == $k } @code_index; my $name = $cinfo{ $k->{key} } ? $cinfo{ $k->{key} }{name} : $k->{crate}; $p .= "* **[$name](./$k->{key}/README.md)**: " . join(', ', map { "[$_->[1]{title}](./$k->{key}/$_->[1]{id}.md)" } @r) . "\n" }
    spit("$root/4-code/README.md", $p);
}
# reports
{
    my $r = "components: $ncomp, lines of Rust in components: $total_loc, classes drawn: $n_classes, public items listed: $n_items\n\ncontainer-level references found in the source (from->to, number of references):\n";
    $r .= "  $_ $crossk{$_}\n" for sort keys %crossk;
    # compare with the container relationships of the DSL
    if ($dsl_file && -f $dsl_file) {
        my $d = slurp($dsl_file); my %model;
        while ($d =~ /^\s*nyx\.(\w+)\s*->\s*nyx\.(\w+)\s+"/mg) { $model{"$1->$2"} = 1 if $cont{$1} && $cont{$2} }
        my @missing = grep { !$model{$_} } sort keys %crossk; my @unseen = grep { !$crossk{$_} } sort keys %model;
        $r .= "\nreferences in the source that the DSL container relationships do not declare:\n" . (@missing ? join('', map { "  $_\n" } @missing) : "  none\n");
        $r .= "\ncontainer relationships of the DSL with no reference in the source (declared only):\n" . (@unseen ? join('', map { "  $_\n" } @unseen) : "  none\n");
    }
    spit("$root/reports/source-scan.txt", $r);
    print $r;
}
# state for phase 2
open my $sf, '>', "$root/reports/.codediagrams.tmp" or die; print $sf join("\0", map { "$_\1$code_diagrams{$_}" } sort keys %code_diagrams); close $sf;
exit 0;

# ---------------------------------------------------------------------------------------------------------- phase 2
sub phase_pages {
    die "phase 2 needs the export folder $root/export/mermaid\n" unless -d "$root/export/mermaid";
    my %code;
    if (-f "$root/reports/.codediagrams.tmp") { my $t = slurp("$root/reports/.codediagrams.tmp"); for my $e (split /\0/, $t) { my ($k, $v) = split /\x01/, $e, 2; $code{$k} = $v } }
    my %mmd; { opendir my $dh, "$root/export/mermaid" or die; for my $f (readdir $dh) { next unless $f =~ /^structurizr-(.+)\.mmd$/; $mmd{$1} = slurp("$root/export/mermaid/$f") } closedir $dh }
    my @viewer;   # [group, title, mermaid]
    my $mm = sub { my ($key) = @_; die "no export for view $key\n" unless $mmd{$key}; return $mmd{$key} };
    my $has_puml = -d "$root/export/plantuml";
    my $dl = sub { my ($key, $up) = @_; return "Export: [Mermaid]($up/export/mermaid/structurizr-$key.mmd)" . ($has_puml ? " · [PlantUML]($up/export/plantuml/structurizr-$key.puml)" : '') };

    # 1 context
    spit("$root/1-context/context.md", "# Level 1: system context\n\nWho uses nyx and what it talks to. Source: the `Context` view of `workspace.dsl`. " . $dl->('Context', '..') . "\n\n```mermaid\n" . $mm->('Context') . "```\n");
    push @viewer, [ '1 Context', 'System context', $mm->('Context') ];
    # 2 containers
    my $t = "# Level 2: containers\n\nOne container per crate of the workspace, plus the files they write. The arrows are exactly the layering table of `crates/nyx-engine/tests/layering.rs`. Source: the `Containers` view. " . $dl->('Containers', '..') . "\n\n```mermaid\n" . $mm->('Containers') . "```\n\n";
    $t .= "| Container | Crate | What it is |\n|---|---|---|\n";
    for my $k (@containers) { my $i = $cinfo{ $k->{key} }; $t .= "| [" . ($i ? $i->{name} : $k->{key}) . "](../3-components/$k->{key}.md) | `crates/$k->{dir}` | " . ($i ? md_esc($i->{desc}) : '') . " |\n" }
    for my $key (sort keys %cinfo) { next if $cont{$key}; $t .= "| " . md_esc($cinfo{$key}{name}) . " | (data) | " . md_esc($cinfo{$key}{desc}) . " |\n" }
    spit("$root/2-containers/containers.md", $t);
    push @viewer, [ '2 Containers', 'Containers and layering', $mm->('Containers') ];
    # 3 components
    for my $k (@containers) {
        my $name = $cinfo{ $k->{key} } ? $cinfo{ $k->{key} }{name} : $k->{crate};
        my $key = "Components_$k->{key}";
        my $p = "# Level 3: components of $name\n\nCrate `crates/$k->{dir}`. " . ($cinfo{ $k->{key} } ? $cinfo{ $k->{key} }{desc} : '') . "\n\nA component is a group of modules with one job; its code page (level 4) has the class diagram and every public item. Arrows to a container are references to that crate's modules (named in the label); arrows between components are `use` and path references in the source. " . $dl->($key, '..') . "\n\n```mermaid\n" . $mm->($key) . "```\n\n";
        $p .= "| Component | Files | Lines | Decisions | What it does |\n|---|---:|---:|---|---|\n";
        for my $c (@{ $k->{comps} }) { my @fl = file_list($k, $c); my $loc = 0; $loc += $_->{lines} for @fl; $p .= "| [$c->{title}](../4-code/$k->{key}/$c->{id}.md) | " . scalar(@fl) . " | $loc | " . md_esc($c->{adr}) . " | " . md_esc($c->{desc}) . " |\n" }
        push @viewer, [ '3 Components', $name, $mm->($key) ];
        for my $v (@{ $k->{views} || [] }) { my $vk = "${key}_$v->{id}"; $p .= "\n## Part: $v->{title}\n\nThe same components, only those of this part (fewer arrows to read). " . $dl->($vk, '..') . "\n\n```mermaid\n" . $mm->($vk) . "```\n"; push @viewer, [ '3 Components', "$name: $v->{title}", $mm->($vk) ] }
        spit("$root/3-components/$k->{key}.md", $p);
    }
    my $idx = "# Level 3: components\n\n" . join('', map { my $n = $cinfo{ $_->{key} } ? $cinfo{ $_->{key} }{name} : $_->{crate}; "* [$n](./$_->{key}.md): " . scalar(@{ $_->{comps} }) . " components\n" } @containers);
    spit("$root/3-components/README.md", $idx);
    # 5 dynamic
    my @dyn = grep { !/^(Context|Containers|Components_)/ } sort keys %mmd;
    my $di = "# Dynamic views\n\nHow the containers work together on one task. Source: the dynamic views of `workspace.dsl`.\n\n";
    for my $key (@dyn) {
        my $desc = $key eq 'SolveFromTheWindow' ? 'What happens when the window solves the active physics.' : $key eq 'AICommand' ? 'An AI host builds a model and solves it through the same commands.' : '';
        spit("$root/5-dynamic/$key.md", "# Dynamic view: $key\n\n$desc " . $dl->($key, '..') . "\n\n```mermaid\n" . $mm->($key) . "```\n");
        $di .= "* [$key]($key.md) $desc\n"; push @viewer, [ '5 Dynamic', $key, $mm->($key) ];
    }
    spit("$root/5-dynamic/README.md", $di);
    for my $k (@containers) { for my $c (@{ $k->{comps} }) { push @viewer, [ "4 Code: " . ($cinfo{ $k->{key} } ? $cinfo{ $k->{key} }{name} : $k->{crate}), $c->{title}, $code{"$k->{key}.$c->{id}"} // "classDiagram\n class missing" ] } }
    # viewer
    my $js = join(",\n", map { my ($g, $ti, $m) = @$_; my $j = $m; $j =~ s/\\/\\\\/g; $j =~ s/"/\\"/g; $j =~ s/\n/\\n/g; $j =~ s/<\//<\\\//g; "[\"" . dq2($g) . "\",\"" . dq2($ti) . "\",\"$j\"]" } @viewer);
    my $html = viewer_html($js, scalar @viewer);
    spit("$root/viewer/index.html", $html);
    print "pages written: " . scalar(@viewer) . " diagrams in the viewer\n";
}
sub dq2 { my $s = shift; $s =~ s/\\/\\\\/g; $s =~ s/"/\\"/g; return $s }
sub viewer_html {
    my ($js, $n) = @_;
    return <<"HTML";
<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>nyx C4 viewer</title>
<style>
:root{--bg:#fff;--fg:#1c2430;--mut:#5b6677;--line:#d5dbe3;--side:#f4f6f9;--acc:#1168bd}
\@media (prefers-color-scheme:dark){:root{--bg:#14181f;--fg:#e6eaf0;--mut:#9aa6b8;--line:#2b3340;--side:#1a2029;--acc:#6cb2ff}}
*{box-sizing:border-box}body{margin:0;font:14px/1.45 system-ui,sans-serif;background:var(--bg);color:var(--fg);display:flex;height:100vh}
nav{width:300px;min-width:240px;overflow:auto;background:var(--side);border-right:1px solid var(--line);padding:12px 8px}
nav h1{font-size:15px;margin:4px 8px 8px}nav h2{font-size:11px;text-transform:uppercase;letter-spacing:.06em;color:var(--mut);margin:14px 8px 4px}
nav a{display:block;padding:3px 8px;border-radius:4px;color:var(--fg);text-decoration:none;cursor:pointer}nav a:hover{background:var(--line)}nav a.on{background:var(--acc);color:#fff}
main{flex:1;overflow:auto;padding:16px 20px}main h2{margin:0 0 2px;font-size:18px}main p{margin:0 0 12px;color:var(--mut)}
#d .edgeLabel,#d .edgeLabel *,#d .label,#d .label *{color:#1c2430 !important;fill:#1c2430}#d .labelBkg,#d .edgeLabel rect{background:#fff !important;fill:#fff !important;opacity:.9}#d svg{max-width:none}#d.fit svg{max-width:100% !important;height:auto}.bar{margin:0 0 8px;color:var(--mut)}
#d{background:#fff;color:#1c2430;border:1px solid var(--line);border-radius:6px;padding:12px;overflow:auto}#err{color:#b3261e;white-space:pre-wrap;font:12px monospace}
\@media (max-width:700px){body{flex-direction:column;height:auto}nav{width:auto;max-height:40vh}}
</style></head><body>
<nav id="nav"><h1>nyx C4 ($n diagrams)</h1></nav>
<main><h2 id="t"></h2><p id="g"></p><div class="bar"><label><input type="checkbox" id="fit" checked> fit to width</label></div><div id="d" class="fit"></div><div id="err"></div></main>
<script src="https://cdn.jsdelivr.net/npm/mermaid\@11.4.1/dist/mermaid.min.js"></script>
<script>
const D=[
$js
];
mermaid.initialize({startOnLoad:false,securityLevel:'loose',maxTextSize:900000,maxEdges:5000,flowchart:{useMaxWidth:false},theme:'default'});
const nav=document.getElementById('nav');let last='';
D.forEach((d,i)=>{if(d[0]!==last){const h=document.createElement('h2');h.textContent=d[0];nav.appendChild(h);last=d[0]}
const a=document.createElement('a');a.textContent=d[1];a.id='a'+i;a.onclick=()=>show(i);nav.appendChild(a)});
async function show(i){document.querySelectorAll('nav a.on').forEach(x=>x.classList.remove('on'));document.getElementById('a'+i).classList.add('on');
document.getElementById('t').textContent=D[i][1];document.getElementById('fit').onchange=e=>document.getElementById('d').classList.toggle('fit',e.target.checked);document.getElementById('g').textContent=D[i][0];document.getElementById('err').textContent='';
try{const r=await mermaid.render('m'+i,D[i][2]);document.getElementById('d').innerHTML=r.svg}catch(e){document.getElementById('d').innerHTML='';document.getElementById('err').textContent=String(e)}location.hash=i}
window.checkAll=async()=>{const bad=[];for(let i=0;i<D.length;i++){try{await mermaid.parse(D[i][2])}catch(e){bad.push([i,D[i][0],D[i][1],String(e).slice(0,300)])}}return bad};
show(Math.max(0,parseInt(location.hash.slice(1))||0));
</script></body></html>
HTML
}
