<?php
$geni = $_GET['geni'];

//controllare se $geni è diverso da vuoto, se è diverso da vuoto non creare file, altrimenti crea file con valori
if(!is_null($geni)){
	//creare file geniSelezionti.txt
	file_put_contents("geniSelezionti.txt",print_r(implode("\n",$geni),true).PHP_EOL );
}

if(empty($_GET['FC']))
	$FC = 0;
	else $FC = $_GET['FC'];

if(empty($_GET['Pvalue']))
	$Pvalue = 0;
	else $Pvalue = $_GET['Pvalue'];

shell_exec("/usr/bin/Rscript Filtri.R $FC $Pvalue 1>result.txt 2>error.log &");

//header("location: result.php");
?>
