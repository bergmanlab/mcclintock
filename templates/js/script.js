function sortTable(id, n, numeric=false){
    console.log(id)
    var table = document.getElementById(id);
    var swapped = true;
    var sortDirection = "lowToHigh";
    var rows,i,j, shouldSwap,x,y,tmp, switchCount = 0;
    rows = table.rows;

    console.log("read table");
    // read in table into an array
    var rowArray = new Array(rows.length-1);
    var valuesArray;
    for (i = 1; i < rows.length; i++){
        valuesArray = new Array(rows[1].getElementsByTagName("td").length+1);
        var valuesLength = rows[1].getElementsByTagName("td").length;
        for (j = 0; j < valuesLength; j++){
            valuesArray[j] = rows[i].getElementsByTagName("td")[j].innerHTML;
        }
        valuesArray[valuesLength] = rows[i].style.display;
        rowArray[i] = valuesArray;
    }

    console.log("sort table");
    // sort array
    while(swapped){
        swapped = false;
        for (i = 1; i < rowArray.length-1; i++){
            x = rowArray[i][n]
            y = rowArray[i+1][n]
            if (sortDirection == "lowToHigh"){
                if (numeric){
                    if (parseInt(x) > parseInt(y)){
                        swapped = true;
                        tmp = rowArray[i];
                        rowArray[i] = rowArray[i+1];
                        rowArray[i+1] = tmp;
                        switchCount += 1;
                        break;
                    }
                } else {
                    if (x.toLowerCase() > y.toLowerCase()){
                        swapped = true;
                        tmp = rowArray[i];
                        rowArray[i] = rowArray[i+1];
                        rowArray[i+1] = tmp;
                        switchCount += 1;
                        break;
                    }
                }

            } else if (sortDirection == "highToLow"){
                if (numeric){
                    if (parseInt(x) < parseInt(y)){
                        swapped = true;
                        tmp = rowArray[i];
                        rowArray[i] = rowArray[i+1];
                        rowArray[i+1] = tmp;
                        switchCount += 1;
                        break;
                    }
                } else {
                    if (x.toLowerCase() < y.toLowerCase()){
                        swapped = true;
                        tmp = rowArray[i];
                        rowArray[i] = rowArray[i+1];
                        rowArray[i+1] = tmp;
                        switchCount += 1;
                        break;
                    }
                }
            }
        }

        if (switchCount == 0 && sortDirection == "lowToHigh"){
            sortDirection = "highToLow";
            swapped = true;
        }
    }

    //clear table
    console.log("clear table");
    for (i = rows.length-1; i > 0; i--){
        table.deleteRow(i);
    }

    console.log("write table");
    // repolulate table with sorted values
    for (i = 1; i < rowArray.length; i++){
        var newRow = table.insertRow(-1);
        var newCell;

        for (j = 0; j < rowArray[1].length-1; j++){
            newCell = newRow.insertCell(j);
            newCell.innerHTML = rowArray[i][j];
        }
        newRow.style.display = rowArray[i][rowArray[1].length-1];
    }

    console.log("fix even-odd column coloring");
    // fix odd even shading
    var rows = table.getElementsByTagName("tr");
    var odd=false;
    for (var i = 1; i < rows.length; i++){
        if (rows[i].style.display == ""){
            if (odd){
                odd = false;
                for (var j = 0; j < rows[i].getElementsByTagName("td").length; j++){
                    rows[i].getElementsByTagName("td")[j].classList.add('even');
                }
            } else {
                odd = true;
                for (var j = 0; j < rows[i].getElementsByTagName("td").length; j++){
                    rows[i].getElementsByTagName("td")[j].classList.remove('even');
                }
            }
        }
    }

}

function filterTableType(inputId, tableId, exactbox, allbox, refbox, nonrefbox){
    var input = document.getElementById(inputId);
    var filter = input.value.toUpperCase();
    var table = document.getElementById(tableId);
    var rows = table.getElementsByTagName("tr");
    var exactMatch = document.getElementById(exactbox).checked;
    var keepAll = document.getElementById(allbox).checked;
    var keepRef = document.getElementById(refbox).checked;
    var keepNonref = document.getElementById(nonrefbox).checked;
    console.log(exactMatch);
    var display = false;
    var data;
    for (var i = 1; i < rows.length; i++){
        display = false;
        data = rows[i].getElementsByTagName("td")[0];
        type = rows[i].getElementsByTagName("td")[1];
        if (data){
            txtValue = data.textContent || data.innerText;
            var typeValue = type.textContent || type.innerText;
            typeValue = typeValue.toUpperCase();

            if (exactMatch){
                if (txtValue.toUpperCase() == filter){
                    if ((typeValue == "ALL" && keepAll) || (typeValue == "REFERENCE" && keepRef) || (typeValue == "NON-REFERENCE" && keepNonref)){
                        display = true;
                    }
                }
            } else {
                if (txtValue.toUpperCase().indexOf(filter) > -1){
                    if ((typeValue == "ALL" && keepAll) || (typeValue == "REFERENCE" && keepRef) || (typeValue == "NON-REFERENCE" && keepNonref)){
                        display = true;
                    }
                }
            }
        }
        if (display){
            rows[i].style.display = "";
        } else {
            rows[i].style.display = "none";
        }
    }
}


function filterTableExtended(input1Id, input2Id, tableId, exactbox, refbox, nonrefbox, filterOptionButton){
    var input1 = document.getElementById(input1Id);
    var filter1 = input1.value.toUpperCase();
    var input2 = document.getElementById(input2Id);
    var filter2 = input2.value.toUpperCase();
    var table = document.getElementById(tableId);
    var rows = table.getElementsByTagName("tr");
    var exactMatch = document.getElementById(exactbox).checked;
    var keepRef = document.getElementById(refbox).checked;
    var keepNonref = document.getElementById(nonrefbox).checked;
    var filterOption = document.getElementById(filterOptionButton).value;
    // console.log(filterOption);
    var match1 = false;
    var match2 = false;
    var data;
    for (var i = 1; i < rows.length; i++){
        match1 = false;
        match2 = false;
        for (var j = 0; j < 2; j++){
            data = rows[i].getElementsByTagName("td")[j];
            type = rows[i].getElementsByTagName("td")[2];
            if (data){
                txtValue = data.textContent || data.innerText;
                var typeValue = type.textContent || type.innerText;
                typeValue = typeValue.toUpperCase();
                if (exactMatch){
                    if (txtValue.toUpperCase() == filter1){
                        if ((typeValue == "REFERENCE" && keepRef) || (typeValue == "NON-REFERENCE" && keepNonref)){
                            match1 = true;
                        }
                    }
                    if (txtValue.toUpperCase() == filter2){
                        if ((typeValue == "REFERENCE" && keepRef) || (typeValue == "NON-REFERENCE" && keepNonref)){
                            match2 = true;
                        }
                    }
                } else {
                    if (txtValue.toUpperCase().indexOf(filter1) > -1){
                        if ((typeValue == "REFERENCE" && keepRef) || (typeValue == "NON-REFERENCE" && keepNonref)){
                            match1 = true;
                        }
                    }
                    if (txtValue.toUpperCase().indexOf(filter2) > -1){
                        if ((typeValue == "REFERENCE" && keepRef) || (typeValue == "NON-REFERENCE" && keepNonref)){
                            match2 = true;
                        }
                    }
                }
            }
        }
        // console.log(match1);
        // console.log(match2);
        if ((filterOption == "or" && (match1 || match2)) || (filterOption == "and" && (match1 && match2))){
            rows[i].style.display = "";
        } else {
            rows[i].style.display = "none";
        }
    }
}

function flipFilterOption(buttonId){
    var button = document.getElementById(buttonId);
    var value = button.value;
    if (value == "and"){
        button.value = "or";
    } else {
        button.value = "and";
    }
}

function hide(divId, hideButtomId){
    var x = document.getElementById(divId);
    var y = document.getElementById(hideButtomId);
    if (x.style.display === "none") {
        x.style.display = "block";
        y.innerHTML = "Hide";
    } else {
        x.style.display = "none";
        y.innerHTML = "Show";
    }
}