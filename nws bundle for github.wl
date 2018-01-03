(* ::Package:: *)

(* ::Input::Initialization:: *)
TDD=MixtureDistribution[{1,1},{UniformDistribution@{-#,#},UniformDistribution@{\[Pi]-#,\[Pi]+#}}]&;


(* ::Input::Initialization:: *)
nanowires=Line[
{{#1,#2},{#3 Cos@#4,#3 Sin@#4}+{#1,#2}}&@@Append[
RandomReal/@{{-5,5},{-5,5},{8,11}},RandomVariate@TDD@(\[Pi]/12)]
]&~Array~50;


(* ::Input::Initialization:: *)
film={{-4.,-10},{4.,10}};
electrodes=Line@{{#1,#2},{#1,-#2}}&@@@film;


(* ::Input::Initialization:: *)
jaggingF=Sort[#~Level~{3}]&/@(
Select[#,RegionWithin[Rectangle@@film,#]&]&/@
Cases[_Point][
l~RegionIntersection~#&/@#1~Join~#2~Complement~{l}
]~Table~{l,#1}~Select~(Length@#>1&)
)~Select~(Length@#>1&)&;


(* ::Input::Initialization:: *)
jagged=jaggingF[nanowires,electrodes];


(* ::Input::Initialization:: *)
grainElectr=Join@@jaggingF[electrodes,nanowires];


(* ::Input::Initialization:: *)
sub=#~AssociationThread~Range@Length@#&[Union@@jagged];
jaglab=MapIndexed[{#1,First@#2}&,jagged/.sub,{2}];


(* ::Input::Initialization:: *)
rsub=KeyMap[{#,_}&,First/@PositionIndex@sub]//Normal;


(* ::Input::Initialization:: *)
{wireedges,grainedges}=
Join@@Map[UndirectedEdge@@#&,Partition[#,2,1]&/@#,{2}]&/@
{jaglab,(Join@@jaglab)~GatherBy~First~Select~(Length@#>1&)};


(* ::Input::Initialization:: *)
activeC=Select[
ConnectedGraphComponents@Graph@Join[wireedges,grainedges],
Length[Flatten[VertexList@#/.rsub]\[Intersection](First/@film)]==2&];


(* ::Input::Initialization:: *)
{wireedges,grainedges}=
#\[Intersection]Join@@EdgeList/@activeC&/@{wireedges,grainedges};


(* ::Input::Initialization:: *)
contactedges=Function[{c,p},
c<->#&/@Select[
VertexList@Graph@Join[wireedges,grainedges],
ContainsAny[#/.rsub,{p},SameTest->Equal]&]
]~MapThread~{{"\!\(\*SuperscriptBox[\(V\), \(+\)]\)","\!\(\*SuperscriptBox[\(V\), \(-\)]\)"},First/@film}//Join@@#&;


(* ::Input::Initialization:: *)
wcontactedges=#->10^6&/@contactedges;


(* ::Input::Initialization:: *)
wgrainedges=Module[{R0=10.,k=1.,n=1./2.},
#->(*1/R0(1+k c^n)*)\[Gamma]g&/@grainedges];


(* ::Input::Initialization:: *)
wwireedges=Module[{R0=10.,k=0.1,n=1.},
#->(*1/R0(1+k c^n)*)\[Gamma]w&/@wireedges];


(* ::Input::Initialization:: *)
\[CapitalOmega][i_,j_,\[ScriptCapitalG]_]:=
WeightedAdjacencyMatrix@\[ScriptCapitalG]//
#-DiagonalMatrix@SparseArray[Total/@#]&//
PseudoInverse//
2#[[i,j]]-#[[i,i]]-#[[j,j]]&


(* ::Input::Initialization:: *)
Graph[Join@@{#1,#2,#3},EdgeStyle->Join[#->{Thick,Dashed,Orange}&/@#1,#->{Black}&/@#2,#->{Thin,Blue}&/@#3](*,VertexLabels\[Rule]"Name"*)]&@@{grainedges,wireedges,contactedges}
