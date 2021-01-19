# Examples

Below are some work though examples to demonstrate how to use the code. The xyz files are the inputs with the POSCAR.vasp files being the below and after, denoted by `_1H2O`.

## CuO4 surface

This is an arbitary 2D material where the Cu was hydrated. This was performed via:
`bin/hydrate_system.x CuO4.xyz 1 10`

## CoO vacany

This is a caculated CoO system (doi: 10.17188/1193863) with 1 Co and 2 O atoms removed. The vacancy was then hydrated with the center denoted by the removed Co's position (labelled as XX). This was performed via:
`bin/hydrate_system.x CoO_vacancy.xyz 60 20`
