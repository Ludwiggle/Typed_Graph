(* ::Package:: *)

(* twin towers distribution function *)
TDD=MixtureDistribution[{1,1},
	{UniformDistribution@{-#,#},UniformDistribution@{\[Pi]-#,\[Pi]+#}}]&;


(* sensing film region *)
film=N@{{-40,-100},{40,100}};
electrodes=Line@{{#1,#2},{#1,-#2}}&@@@film;


(* jagged array of intersection function *)
jaggingF=Sort[#~Level~{3}]&/@(
	Select[#,RegionWithin[Rectangle@@film,#]&]&/@
	Cases[_Point][
	l~RegionIntersection~#&/@#1~Join~#2~Complement~{l}
	]~Table~{l,#1}~Select~(Length@#>1&)
	)~Select~(Length@#>1&)&;


(* neck-controlled grain-boundary-controlled chemiresitance *)
rN[Nr_,L_,aN_,aB_]:=(1/(aN ((-1+L/2+Nr)^2-E^(-aB (-1+Nr)^2) (-1+Nr) (-1+L+Nr))))
rGB[Nr_,L_,a_,b_]:=(b E^(a (-1+Nr)^2))/L^2


(* resistance distance *)
om[i_,j_,G_]:=
WeightedAdjacencyMatrix@G//
#-DiagonalMatrix@SparseArray[Total/@#]&//
PseudoInverse//
2#[[i,j]]-#[[i,i]]-#[[j,j]]&

(* number of nanowires *)
nn=20;
(* nn random nanowires *)
nanowires=Line[
	{{#1,#2},{#3 Cos@#4,#3 Sin@#4}+{#1,#2}}&@@Append[
	RandomReal/@{{-50,50},{-50,50},{80,110}},RandomVariate@TDD@(\[Pi]/2)]
	]&~Array~nn;


(* jagged arrays of intersections *)
jagged=jaggingF[nanowires,electrodes];
grainElectr=Join@@jaggingF[electrodes,nanowires];


(* substitution rule {x,y}\[Rule]{gbID,nw}  *)
sub=#~AssociationThread~Range@Length@#&[Union@@jagged];
(* reverse substitution {nw, gbID}\[Rule]{x,y} *)
rsub=KeyMap[{#,_}&,First/@PositionIndex@sub]//Normal;


(* jagged array of vertex vectors *)
jaglab=MapIndexed[{#1,First@#2}&,jagged/.sub,{2}];


(* edges *)
{wireedges,grainedges}=Function[edges,
	#\[Intersection]Join@@EdgeList/@Select[
	ConnectedGraphComponents@Graph@(Join@@edges),
	Length[Flatten@Rationalize[VertexList@#/.rsub,10^(-6)]\[Intersection]Rationalize[First/@film,10^(-6)]]==2&]&/@edges][
	Join@@Map[UndirectedEdge@@#&,Partition[#,2,1]&/@#,{2}]&/@
	{jaglab,(Join@@jaglab)~GatherBy~First~Select~(Length@#>1&)}];

contactedges=Function[{c,p},
	c\[UndirectedEdge]#&/@Select[
	VertexList@Graph@Join[wireedges,grainedges],
	ContainsAny[#/.rsub,{p},SameTest->Equal]&]
	]~MapThread~{{"V^+","V^-"},First/@film}//Join@@#&;



(* weighted edges *)
wcontactedges=#->10^4&/@contactedges;

wwireedges[c_,diam_]:=
	#1->1/(#2 rN[c,diam,1.069,4.39])&~MapThread~{
	wireedges,Norm@@@Differences/@List@@@(wireedges/.rsub)};

wgrainedges[c_,diam_]:=#->1/rGB[c,diam,3.66,5]&/@grainedges;


(* final unweighted graph *)
finalg=Graph[
	Join@@{#1,#2,#3},
	EdgeStyle->Join[
		#->{Thick,Dashed,Orange}&/@#1,
		#->{Black}&/@#2,
		#->{Thin,Blue}&/@#3]]&@@{grainedges,wireedges,contactedges};


(* weighted graph (Typed Graph) function *)
tg[c_,diam_]:=Graph[finalg,EdgeWeight->Join[
	wgrainedges[c,diam],wwireedges[c,diam],wcontactedges]];


(* response to c=0.5 gas surface coverage
 for diameter d=7.7 times the depletion depth *)
respose=Module[{c=.5,d=7.7,tg0,tg1,r0,r1,prb},
	{tg0,tg1}=tg[#,d]&/@{0.,c};
	prb=Position[VertexList@tg0,#]&/@{"V^+","V^-"}//Flatten;
	r0=om[prb[[1]],prb[[2]],tg0];
	r1=om[prb[[1]],prb[[2]],tg1];
	r0/r1]~Check~"no short circuit"//Quiet
