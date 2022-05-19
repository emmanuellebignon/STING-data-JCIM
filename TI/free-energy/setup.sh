#!/bin/sh
#
# Setup for the free energy simulations: creates and links to the input file as
# necessary.  Two alternative for the de- and recharging step can be used.
#


basedir=../setup
top=$(pwd)
setup_dir=$(cd "$basedir"; pwd)

for system in protein complex; do
  if [ \! -d $system ]; then
    mkdir $system
  fi

  cd $system

  for w in 0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0; do
    if [ \! -d $w ]; then
      mkdir $w
    fi

    sed -e "s/%L%/$w/" $top/min.tmpl > $w/min.in
    sed -e "s/%L%/$w/" $top/heat.tmpl > $w/heat.in
    sed -e "s/%L%/$w/" $top/prod.tmpl > $w/ti.in

    (
      cd $w
      ln -sf $setup_dir/merged_${system}.parm7 ti.parm7
      ln -sf $setup_dir/merged_${system}.rst7  ti.rst7
    )
  done

  cd $top
done

