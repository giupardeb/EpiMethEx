<html>
  <script
    src="https://code.jquery.com/jquery-2.2.4.min.js"
    integrity="sha256-BbhdlvQf/xTY9gja0Dq3HiwQF8LaCRTXxZKRutelT44="
    crossorigin="anonymous"></script>

  <?php

  $obj = json_decode(file_get_contents("UPvsMIDordinato.json"), true);
  $string = "";
  $duplicatecg = array();

  echo '<table border="1">
      <tr>
          <td colspan="4" style="text-align: center;"><h3>UPvsMID</h3></td>
        </tr>
        <tr>
          <td style="text-align: center;"><b>Gene</b></td>
          <td style="text-align: center;"><b>Fold Change</b></td>
          <td style="text-align: center;"><b>P Value</b></td>
          <td style="text-align: center;"><b>CG</b></td>
        </tr>';

  foreach ( $obj as $k=>$v ) {   
    foreach ( $obj[$k] as $z=>$j ) {
      echo '<tr><td>'.$obj[$k][$z]['Gene'].'</td>';
      echo '<td>'.$obj[$k][$z]['FC'].'</td>';
      echo '<td>'.$obj[$k][$z]['pValue'].'</td>';

      $posizione = explode(';;', $obj[$k][$z]['Posizione']);
      
      echo '<td>&nbsp;<a href="#'.$obj[$k][$z]['Gene'].'" onClick="show_hide(this);" id="'.$obj[$k][$z]['Gene'].'">Show them all</a>&nbsp;';
      echo '<p><div id="'.$obj[$k][$z]['Gene'].'_ul" style="display: none;"><ul>';
      for($i = 0; $i < count($posizione); $i = $i+3) {
    
        if($string=="")
          $string = $posizione[$i].'('.$posizione[$i+1].','.$posizione[$i+2];
        else
          $string = $string.','.$posizione[$i+1].','.$posizione[$i+2];
      
        if($posizione[$i] != $posizione[$i+3]) {
          echo '<li>'.$string.') </li>';
          $string = "";
        } else { //visualizzo all'utente i CG duplicati
          if(!in_array($posizione[$i], $duplicatecg, true))
            array_push($duplicatecg, $posizione[$i]);
        }
      }
      echo '</ul></div></p></td>';

    }
  }

  /*echo '<tfoot>
      <tr>
        <td colspan="2">CG duplicati <pre>'; print_r($duplicatecg); echo '</pre></td>
      </tr>
    <tfoot>';*/
  echo '</table>';

  ?>

  <script>
  function show_hide(param) {
    var $id = $(param).attr("id")+"_ul";
    $("#"+$id).toggle();
  }
  </script>
</html>