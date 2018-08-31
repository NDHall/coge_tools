**__**notes on qac files validating the results.

* Connection network taken from syntenic block links take from dag chainer output and combined with bedtools overlap
here we are validating the output. One thing I noticed that we were getting connections more than 4 with our 2:1 
quota align clustering. This is not a fault of quota align. I am using a differnt approach to look for over lap than 
the extremely capable authors of quota align. 

* histogram of length of connected nodes, this should really be between 2-4. In the best of all possible worlds 3. 

![alt text][histogram] 

![alt text][subgraph1]

```
sub-graph contains:

'Super-Scaffold_69_567283_835551': 'scaffold197_295422_570637', 
'scaffold133_434808_490887': 'Super-Scaffold_307_2523554_2587946',
'scaffold133_428029_492266': 'Super-Scaffold_2499_1894054_1977853',
'scaffold197_294964_572620': 'Super-Scaffold_2499_1644284_1903681'
 ```
* check parsing from initial file. 
```bash
egrep -w 'Super-Scaffold_69' 51576_52024.CDS-CDS.last.tdd10.cs0.filtered.dag.all.go_D20_g10_A5.aligncoords.qac2.1.40  |\
egrep -w 'scaffold197'| cut -d'|' -f -6
a51576_Super-Scaffold_69	Super-Scaffold_69||567284||568612| #starting coordinate correct
a51576_Super-Scaffold_69	Super-Scaffold_69||584051||587263|
a51576_Super-Scaffold_69	Super-Scaffold_69||591789||594845|
a51576_Super-Scaffold_69	Super-Scaffold_69||598428||603031|
a51576_Super-Scaffold_69	Super-Scaffold_69||605616||607908|
a51576_Super-Scaffold_69	Super-Scaffold_69||617923||620011|
a51576_Super-Scaffold_69	Super-Scaffold_69||620965||622244|
a51576_Super-Scaffold_69	Super-Scaffold_69||626254||637617|
a51576_Super-Scaffold_69	Super-Scaffold_69||639467||639655|
a51576_Super-Scaffold_69	Super-Scaffold_69||646938||650463|
a51576_Super-Scaffold_69	Super-Scaffold_69||651681||657557|
a51576_Super-Scaffold_69	Super-Scaffold_69||659464||663690|
a51576_Super-Scaffold_69	Super-Scaffold_69||668275||670523|
a51576_Super-Scaffold_69	Super-Scaffold_69||678263||680528|
a51576_Super-Scaffold_69	Super-Scaffold_69||685080||689134|
a51576_Super-Scaffold_69	Super-Scaffold_69||692144||693256|
a51576_Super-Scaffold_69	Super-Scaffold_69||697127||698011|
a51576_Super-Scaffold_69	Super-Scaffold_69||726270||728463|
a51576_Super-Scaffold_69	Super-Scaffold_69||736320||739222|
a51576_Super-Scaffold_69	Super-Scaffold_69||739800||742199|
a51576_Super-Scaffold_69	Super-Scaffold_69||743283||746246|
a51576_Super-Scaffold_69	Super-Scaffold_69||747152||751386|
a51576_Super-Scaffold_69	Super-Scaffold_69||767029||768340|
a51576_Super-Scaffold_69	Super-Scaffold_69||768872||774698|
a51576_Super-Scaffold_69	Super-Scaffold_69||776612||777942|
a51576_Super-Scaffold_69	Super-Scaffold_69||779116||781871|
a51576_Super-Scaffold_69	Super-Scaffold_69||782680||785582|
a51576_Super-Scaffold_69	Super-Scaffold_69||786055||787917|
a51576_Super-Scaffold_69	Super-Scaffold_69||804038||805087|
a51576_Super-Scaffold_69	Super-Scaffold_69||810409||812803|
a51576_Super-Scaffold_69	Super-Scaffold_69||821935||823519|
a51576_Super-Scaffold_69	Super-Scaffold_69||823617||825724|
a51576_Super-Scaffold_69	Super-Scaffold_69||835093||835551| #ending coordinate correct

# second half
egrep -w 'Super-Scaffold_69' 51576_52024.CDS-CDS.last.tdd10.cs0.filtered.dag.all.go_D20_g10_A5.aligncoords.qac2.1.40 |\
egrep -w 'scaffold197'| awk '{print $5,$6,$7}' |cut -d'|' -f -6
b52024_scaffold197 scaffold197||570637||571227| # notice this one is in the opposite direction
b52024_scaffold197 scaffold197||561999||564885|
b52024_scaffold197 scaffold197||553127||555313|
b52024_scaffold197 scaffold197||544783||549628|
b52024_scaffold197 scaffold197||541294||543496|
b52024_scaffold197 scaffold197||527276||529222|
b52024_scaffold197 scaffold197||524742||526050|
b52024_scaffold197 scaffold197||510218||521238|
b52024_scaffold197 scaffold197||508168||508611|
b52024_scaffold197 scaffold197||499410||503434|
b52024_scaffold197 scaffold197||491114||496486|
b52024_scaffold197 scaffold197||477248||488829|
b52024_scaffold197 scaffold197||445079||447562|
b52024_scaffold197 scaffold197||434750||438752|
b52024_scaffold197 scaffold197||426767||431480|
b52024_scaffold197 scaffold197||417802||418911|
b52024_scaffold197 scaffold197||412635||414192|
b52024_scaffold197 scaffold197||405313||408627|
b52024_scaffold197 scaffold197||390946||396530|
b52024_scaffold197 scaffold197||387745||390402|
b52024_scaffold197 scaffold197||383793||385549|
b52024_scaffold197 scaffold197||379156||382887|
b52024_scaffold197 scaffold197||362724||364033|
b52024_scaffold197 scaffold197||356696||362217|
b52024_scaffold197 scaffold197||341761||343265|
b52024_scaffold197 scaffold197||337663||339105|
b52024_scaffold197 scaffold197||333879||336828|
b52024_scaffold197 scaffold197||330993||332870|
b52024_scaffold197 scaffold197||317197||317478|
b52024_scaffold197 scaffold197||308865||312549|
b52024_scaffold197 scaffold197||302157||303029|
b52024_scaffold197 scaffold197||299913||302058|
b52024_scaffold197 scaffold197||294965||295423| # no problem parser still takes 1 from starting coordinate to make it 
                                                # bed coordinates.
```
 * The next connection is between the parent genomes. This should either occur through identical bed coordinates, or more
   likely through the overlap of the regions. Here we see parent is connected to parent through overlap. we won't keep 
   printing entire blocks of code. but the overlap is correct here too. 
 ```bash
 egrep -w 'Super-Scaffold_2499' 51576_52024.CDS-CDS.last.tdd10.cs0.filtered.dag.all.go_D20_g10_A5.aligncoords.qac2.1.40\
 | egrep -w 'scaffold197'| cut -d'|' -f -6
a51576_Super-Scaffold_2499	Super-Scaffold_2499||1644285||1644743| #coordinates are correctly adjusted again.
...
a51576_Super-Scaffold_2499	Super-Scaffold_2499||1901141||1903681|
```  

* We can continue to trust the underlying parsing. We can now trace overlaps based on coordinate overlaps. 
    for example
    Super-Scaffold_2499_1894054_1977853': 'Super-Scaffold_2499_1644284_1903681' 
    we can take the overlap as follows
    ```
    bash
    echo $((1977853 - 1903681)) 
    74172
    ```
  an overlap ~70K could easily lead us to measure same region twice. 
  ```
  bash
    >egrep -w 'Super-Scaffold_2499' \
    51576_52024.CDS-CDS.last.tdd10.cs0.filtered.dag.all.go_D20_g10_A5.aligncoords.qac2.1.40 |\
    >egrep -w 'scaffold197'| cut -d'|' -f -8  >tmp 
    # now get second region 
    >egrep -w 'Super-Scaffold_2499' \
    51576_52024.CDS-CDS.last.tdd10.cs0.filtered.dag.all.go_D20_g10_A5.aligncoords.qac2.1.40 |\
    >egrep -w 'scaffold133'| cut -d'|' -f -8>> tmp # notice this is a different parent genome
    cut -d'|' -f 7- tmp | sort| uniq -c 
      1 PR202_augustus_masked-Super-Scaffold_2499-processed-gene-1.430-mRNA-1|
      1 PR202_augustus_masked-Super-Scaffold_2499-processed-gene-1.433-mRNA-1|
      1 PR202_augustus_masked-Super-Scaffold_2499-processed-gene-1.434-mRNA-1|
      1 PR202_augustus_masked-Super-Scaffold_2499-processed-gene-1.435-mRNA-1|
      1 PR202_augustus_masked-Super-Scaffold_2499-processed-gene-1.436-mRNA-1|
      1 PR202_augustus_masked-Super-Scaffold_2499-processed-gene-1.442-mRNA-1|
      1 PR202_augustus_masked-Super-Scaffold_2499-processed-gene-1.443-mRNA-1|
      1 PR202_augustus_masked-Super-Scaffold_2499-processed-gene-1.446-mRNA-1|
      1 PR202_augustus_masked-Super-Scaffold_2499-processed-gene-1.454-mRNA-1|
      1 PR202_augustus_masked-Super-Scaffold_2499-processed-gene-1.456-mRNA-1|
      1 PR202_augustus_masked-Super-Scaffold_2499-processed-gene-1.461-mRNA-1|
      1 PR202_augustus_masked-Super-Scaffold_2499-processed-gene-1.463-mRNA-1|
      1 PR202_augustus_masked-Super-Scaffold_2499-processed-gene-1.548-mRNA-1|
      1 PR202_augustus_masked-Super-Scaffold_2499-processed-gene-1.558-mRNA-1|
      1 PR202_augustus_masked-Super-Scaffold_2499-processed-gene-1.563-mRNA-1|
      1 PR202_augustus_masked-Super-Scaffold_2499-processed-gene-1.564-mRNA-1|
      1 PR202_augustus_masked-Super-Scaffold_2499-processed-gene-1.565-mRNA-1|
      2 PR202_augustus_masked-Super-Scaffold_2499-processed-gene-1.572-mRNA-1| # this one
      1 PR202_maker-Super-Scaffold_2499-augustus-gene-1.632-mRNA-1|
      1 PR202_maker-Super-Scaffold_2499-augustus-gene-1.635-mRNA-1|
      1 PR202_maker-Super-Scaffold_2499-augustus-gene-1.636-mRNA-1|
      1 PR202_maker-Super-Scaffold_2499-augustus-gene-1.643-mRNA-1|
      1 PR202_maker-Super-Scaffold_2499-augustus-gene-1.644-mRNA-1|
      1 PR202_maker-Super-Scaffold_2499-augustus-gene-1.647-mRNA-1|
      1 PR202_maker-Super-Scaffold_2499-augustus-gene-1.648-mRNA-1|
      1 PR202_maker-Super-Scaffold_2499-augustus-gene-1.651-mRNA-1|
      1 PR202_maker-Super-Scaffold_2499-augustus-gene-1.692-mRNA-1|
      1 PR202_maker-Super-Scaffold_2499-augustus-gene-1.693-mRNA-1|
      1 PR202_maker-Super-Scaffold_2499-augustus-gene-1.695-mRNA-1|
      1 PR202_maker-Super-Scaffold_2499-augustus-gene-1.696-mRNA-1|
      1 PR202_maker-Super-Scaffold_2499-augustus-gene-1.699-mRNA-1|
      1 PR202_maker-Super-Scaffold_2499-augustus-gene-1.700-mRNA-1|
      1 PR202_maker-Super-Scaffold_2499-augustus-gene-1.704-mRNA-1|
      1 PR202_maker-Super-Scaffold_2499-augustus-gene-1.705-mRNA-1|
      1 PR202_maker-Super-Scaffold_2499-augustus-gene-1.706-mRNA-1|
      1 PR202_maker-Super-Scaffold_2499-augustus-gene-1.709-mRNA-1|
      1 PR202_maker-Super-Scaffold_2499-augustus-gene-1.712-mRNA-1|
    # sure enough 
    >egrep -w 'PR202_augustus_masked-Super-Scaffold_2499-processed-gene-1.572-mRNA-1' \
    51576_52024.CDS-CDS.last.tdd10.cs0.filtered.dag.all.go_D20_g10_A5.aligncoords.qac2.1.40 | awk '{print $1,$5}' 

        a51576_Super-Scaffold_2499 b52024_scaffold197
        a51576_Super-Scaffold_2499 b52024_scaffold133
    ```
  
* a couple more interesting examples with node labels.

![alt text][subgraph2]
```
'scaffold202_158037_251120', 'Super-Scaffold_140_817200_947562', 
'Super-Scaffold_88_464670_702466', 'scaffold241_198212_322082', 
'scaffold120_83390_244358', 'Super-Scaffold_2628_1903023_2128530', 
'Super-Scaffold_2628_1989132_2439739', 'Super-Scaffold_2628_1731198_1933381', 
'scaffold202_181873_251120', 'Super-Scaffold_88_4619_268690', 
'scaffold241_195908_327789', 'Super-Scaffold_39_5095216_5389077', 
'Super-Scaffold_95_1246594_1564111', 'scaffold120_71452_176589', 
'scaffold120_152824_512425'}
```

*one last example I liked how it looked. 

![alt text][subgraph3]
```bash
  'Super-Scaffold_134_315978_632853', 'Super-Scaffold_134_655195_795329', 
  'Super-Scaffold_146_790348_944120', 'scaffold666_79114_170649', 
  'Super-Scaffold_134_288383_473856', 'scaffold47_19879_138493',
  'Super-Scaffold_128_2952039_3142249', 'scaffold13_348645_841629',
  'scaffold467_17619_253548', 'scaffold47_21966_136916', 
  'Super-Scaffold_128_3106960_3508023', 'Super-Scaffold_128_2668821_2783585', 
  'Super-Scaffold_152_8522366_8678424', 'scaffold91_109862_397501', 
  'scaffold91_109539_401161', 'scaffold1449_39340_150403', 
  'scaffold91_53216_124856', 'scaffold666_78084_176553', 
  'scaffold13_439383_849553', 'Super-Scaffold_146_670357_838725', 
  'Super-Scaffold_134_472380_696152', 'Super-Scaffold_311_3459784_3631500', 
  'Super-Scaffold_128_2794179_3115952', 'Super-Scaffold_134_9741_322518', 
  'scaffold1449_11629_150403'

```
[histogram]:hist_figure_1.png "histogram of connection per sub-graph"
[subgraph1]:Figure_1.png "example 1"
[subgraph2]:Figure_2.png "example 2"
[subgraph3]:Figure_3.png "example 3"