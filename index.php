<html>
<head>
<title> PHP Test Script </title>
</head>
<body>

<?php
echo("<form action='/script.php'>");
echo("<select name='geni[]' multiple size='20'>");

	$handle = fopen("geni.txt", "r");
	if ($handle) {
    	while (($line = fgets($handle)) !== false) {
         echo("<option value=".$line.">".$line."</option>");
    	}

    fclose($handle);
	} else {
    // error opening the file.
	}
echo("</select>");
echo("<br><br>");
echo("<input type='number' step=any name='FC' placeholder='FoldChange'><br><br>");
echo("<input type='number' step=any name='Pvalue' placeholder='P_value'><br><br>");


echo("<input type='submit'>");
echo("</form>");

?>

</body>
</html>
