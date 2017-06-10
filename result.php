<?php

//echo "<meta http_equiv='refresh' content='5'";
$string = file_get_contents("test.json");
$obj = json_decode($string, true);

$flag = false;
echo "<pre>" . $string . "</pre>";

/*
foreach ( $obj as $k=>$v ) {
  
  echo $obj[$k]['name'].'<br />';
  
  if( is_array($obj[$k]['CG']) == 1 ) { 
	foreach ( $obj[$k]['CG'] as $j=>$z ) {
		foreach($z as $x=>$y) {
			if($flag == false) {
				echo $j.' :'.$y.' '; 
				$flag = true;
			}
			else {
				echo $y.' ';
			}
		}
	 echo '<br />';
     $flag = false;
  	}
  } 
  else echo "NOOO ARRAY"."<br />";
}
*/

?>
