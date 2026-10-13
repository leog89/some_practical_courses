libname mybib "/home/u64137191/UE25";

/* Importation des jeux de données à étudier */

PROC IMPORT OUT=mybib.molecule 
	DATAFILE="/home/u64137191/UE25/molecule.csv"
	DBMS=CSV;
	GETNAMES=YES;
RUN;

PROC CONTENTS DATA=mybib.molecule;
RUN;

PROC IMPORT OUT=mybib.ozone
	DATAFILE="/home/u64137191/UE25/ozone.xls"
	DBMS=XLS;
	GETNAMES=YES;
RUN;

PROC CONTENTS DATA=mybib.ozone;
RUN;

PROC IMPORT OUT=mybib.cancerp
	DATAFILE="/home/u64137191/UE25/cancerp.xls"
	DBMS=XLS;
	GETNAMES=YES;
RUN;

PROC CONTENTS DATA=mybib.cancerp;
RUN;

/* 					Exercice 1 					*/

/* Statistiques descriptives */

PROC MEANS DATA=mybib.molecule NOLABELS;
	VAR mesure;
	TITLE "Statistiques descriptives";
	CLASS traitement;
RUN;

/* Test de normalité */

PROC UNIVARIATE DATA=mybib.molecule NORMAL;
	VAR mesure;
	BY traitement;
	QQPLOT mesure / NORMAL (MU=EST SIGMA=EST);
RUN;

/* Test de Fisher + T-test unilatéral */

PROC TTEST DATA=mybib.molecule SIDE = L;
	VAR mesure;
	CLASS traitement;
RUN;


/* 				 Exercice 2						 */

PROC MEANS DATA=mybib.ozone NOLABELS;
	VAR maxO3 T9 T12 T15 Ne9 Ne12 Ne15 Vx9 Vx12 Vx15;
RUN;

/* Centrage-réduction */

PROC STANDARD DATA=mybib.ozone MEAN=0 STD=1 out=mybib.ozone;
    VAR T9 Ne9 Vx9 T12 Ne12 Vx12 T15 Ne15 Vx15 ;
RUN;

/* Vérification de la relation entre Y (ozone) et chaque var explicative */

PROC SGPLOT DATA=mybib.ozone;
	scatter x=T9 y=maxO3;
RUN;

PROC SGPLOT DATA=mybib.ozone;
	scatter x=T12 y=maxO3;
RUN;

PROC SGPLOT DATA=mybib.ozone;
	scatter x=T15 y=maxO3;
RUN;

PROC SGPLOT DATA=mybib.ozone;
	scatter x=Ne9 y=maxO3;
RUN;

PROC SGPLOT DATA=mybib.ozone;
	scatter x=Ne12 y=maxO3;
RUN;

PROC SGPLOT DATA=mybib.ozone;
	scatter x=Ne15 y=maxO3;
RUN;

PROC SGPLOT DATA=mybib.ozone;
	scatter x=Vx9 y=maxO3;
RUN;

PROC SGPLOT DATA=mybib.ozone;
	scatter x=Vx12 y=maxO3;
RUN;

PROC SGPLOT DATA=mybib.ozone;
	scatter x=Vx15 y=maxO3;
RUN;

/* résidus etc */

PROC REG DATA=mybib.ozone;
	model maxO3 = T9 T12 T15 Ne9 Ne12 Ne15 Vx9 Vx12 Vx15;
RUN;

PROC UNIVARIATE DATA=mybib.ozone NORMAL;
	VAR resreg;
RUN;

/* selection par vif < 3 */

PROC REG DATA=mybib.ozone;
	model maxO3 = T9 T12 T15 Ne9 Ne12 Ne15 Vx9 Vx12 Vx15 / VIF;
RUN;

PROC REG DATA=mybib.ozone;
	model maxO3 = T9 T15 Ne9 Ne12 Ne15 Vx9 Vx12 Vx15 / VIF;
RUN;

PROC REG DATA=mybib.ozone;
	model maxO3 = T9 Ne9 Ne12 Ne15 Vx9 Vx12 Vx15 / VIF;
RUN;

PROC REG DATA=mybib.ozone;
	model maxO3 = T9 Ne9 Ne12 Ne15 Vx9 Vx15 / VIF;
RUN;

PROC REG DATA=mybib.ozone;
	model maxO3 = T9 Ne9 Ne15 Vx9 Vx15 / VIF;
RUN;

/* Selection backward sur la base de la significativité des coeffs */

PROC REG DATA=mybib.ozone;
	model maxO3 = T9 Ne9 Ne15 Vx9 Vx15 / selection = backward slstay = 0.05;
RUN;

/* Modèle final*/

PROC REG DATA=mybib.ozone;
	model maxO3 = T9 Ne9 Vx9;
	OUTPUT OUT = mybib.ozoneres r = resreg;
RUN;

PROC REG DATA=mybib.ozone plots=diagnostics(unpack);
	model maxO3 = T9 Ne9 Vx9;
RUN;

/*					Exercice 3					*/

/* Centrage-réduction des variables quantitatives */

PROC STANDARD DATA=mybib.cancerp MEAN=0 STD=1 out=mybib.cancerp;
    VAR age acide logacid;
RUN;


PROC GENMOD DATA=mybib.cancerp; /* DESCENDING pour modéliser Y=1 et non Y=0 */
	MODEL Y = age acide rayonx taille grade logacid / DIST=binomial LINK=logit;
RUN;

/* Processus de sélection du modèle optimal */

/* Sélection sur la base de la multicolinéarité entre les variables avec le VIF > 3 */
/* Pour pouvoir utiliser les VIF, on doit se ramener à une RLM */

PROC REG DATA=mybib.cancerp;
	MODEL Y = age acide rayonx taille grade logacid / VIF;
RUN;

PROC REG DATA=mybib.cancerp;
	MODEL Y = age acide rayonx taille grade / VIF;
RUN;

/* 2) Sélection backward par LLRT */

PROC GENMOD DATA=mybib.cancerp DESCENDING;
	MODEL Y = age acide rayonx taille grade / DIST=binomial LINK=logit type3;
RUN;

PROC GENMOD DATA=mybib.cancerp DESCENDING; 
	MODEL Y = age acide rayonx taille / DIST=binomial LINK=logit type3;
RUN;

PROC GENMOD DATA=mybib.cancerp DESCENDING;
	MODEL Y = acide rayonx taille / DIST=binomial LINK=logit type3;
RUN;

PROC GENMOD DATA=mybib.cancerp DESCENDING; 
	MODEL Y = rayonx taille / DIST=binomial LINK=logit type3;
RUN;

/* Modèle optimal sans interaction */

PROC GENMOD DATA=mybib.cancerp DESCENDING; 
	MODEL Y = rayonx taille / DIST=binomial LINK=logit;
RUN;

/* Test si l'interaction est significative */

PROC GENMOD DATA=mybib.cancerp DESCENDING; 
	MODEL Y = rayonx taille rayonx*taille / DIST=binomial LINK=logit type3;
RUN;

PROC GENMOD DATA=mybib.cancerp DESCENDING; 
	MODEL Y = rayonx taille / DIST=binomial LINK=logit type3;
RUN;