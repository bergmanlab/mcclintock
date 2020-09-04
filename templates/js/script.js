function sortTable(id, n, numeric=false){
    var t0 = performance.now();
    console.log(id)
    var table = document.getElementById(id);
    var rows, i, j;
    rows = table.rows;

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

    var t1 = performance.now();
    console.log("read table: " + (t1-t0) + "ms");

    invert = isSorted(rowArray, n, numeric=numeric);
    if (!invert){
        invert = isSorted(rowArray, n, numeric=numeric, inverted=true);
    }
    var t2 = performance.now();
    console.log("check if already sorted: "+ (t2-t1) + "ms");

    if(invert){
        rowArray = invertArray(rowArray);
        var t3 = performance.now();
        console.log("invert: " + (t3-t2) + "ms");
    } else{

        sortByColumn(rowArray, n, numeric=numeric);
        var t3 = performance.now();
        console.log("sort: " + (t3-t2) + "ms");
    }


    //clear table
    table.style.display = "none";
    for (i = rows.length-1; i > 0; i--){
        table.deleteRow(i);
    }

    var t4 = performance.now();
    console.log("clear table: " + (t4-t3) + "ms");

    // repolulate table with sorted values
    for (i = 1; i < rowArray.length; i++){
        var newRow = document.createElement("tr");
        var newCell;

        for (j = 0; j < rowArray[1].length-1; j++){
            newCell = newRow.insertCell(j);
            newCell.innerHTML = rowArray[i][j];
            newCell.classList.add("values");
        }
        newRow.style.display = rowArray[i][rowArray[1].length-1];

        table.tBodies[0].appendChild(newRow);
    }

    var t5 = performance.now();
    console.log("write table: " + (t5-t4) + "ms");

    table.style.display = "";
    var t7 = performance.now();
    console.log("total duration: "+ (t7-t0) + "ms");
}

function isSorted(arr, col, numeric=false, inverted=false){
    for (var i =1; i < arr.length-1; i++){
        if (numeric){
            if (inverted){
                if (parseInt(arr[i][col]) < parseInt(arr[i+1][col])){
                    return false;
                }
            } else {
                if (parseInt(arr[i][col]) > parseInt(arr[i+1][col])){
                    return false;
                }
            }

        } else {
            if (inverted){
                if (arr[i][col].toLowerCase() < arr[i+1][col].toLowerCase()){
                    return false;
                }
            } else {
                if (arr[i][col].toLowerCase() > arr[i+1][col].toLowerCase()){
                    return false;
                }
            }

        }
    }
    return true;
}

function invertArray(arr){
    var invertedArray = [...arr];
    for (var i = 1; i < arr.length; i++){
        invertedArray[i] = arr[(arr.length)-i]
    }

    return invertedArray;
}

function sortByColumn(arr, colIdx, numeric=false){
    if (numeric){
        function sortFunction(a,b){
            if (parseInt(a[colIdx]) > parseInt(b[colIdx])){
                return 1;
            } else {
                return -1;
            }
        }
    } else {
        function sortFunction(a,b){
            if (a[colIdx].toLowerCase() > b[colIdx].toLowerCase()){
                return 1;
            } else {
                return -1;
            }
        }
    }

    arr.sort(sortFunction);
    tmp = arr[-1];
    arr.unshift(tmp);
    arr.pop();
    return arr;
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
                if ((txtValue.toUpperCase() == filter) || (filter === "")){
                    if ((typeValue == "ALL" && keepAll) || (typeValue == "REFERENCE" && keepRef) || (typeValue == "NON-REFERENCE" && keepNonref)){
                        display = true;
                    }
                }
            } else {
                if ((txtValue.toUpperCase().indexOf(filter) > -1)  || (filter === "")){
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

    table.style.display = "none";
    var t0 = performance.now();

    var rowArray = new Array(rows.length);
    var valuesArray;
    for (var i = 1; i < rows.length; i++){
        match1 = false;
        match2 = false;

        valuesArray = new Array(rows[1].getElementsByTagName("td").length+1);

        for (var j = 0; j < rows[1].getElementsByTagName("td").length; j++){
            valuesArray[j] = rows[i].getElementsByTagName("td")[j].innerHTML;
        }

        for (var j = 0; j < 2; j++){
            data = rows[i].getElementsByTagName("td")[j];
            type = rows[i].getElementsByTagName("td")[2];
            if (data){
                txtValue = data.textContent || data.innerText;
                var typeValue = type.textContent || type.innerText;
                typeValue = typeValue.toUpperCase();
                if (exactMatch){
                    if (txtValue.toUpperCase() == filter1 || filter1 == ""){
                        if ((typeValue == "REFERENCE" && keepRef) || (typeValue == "NON-REFERENCE" && keepNonref)){
                            match1 = true;
                        }
                    }
                    if (txtValue.toUpperCase() == filter2 || filter2 == ""){
                        if ((typeValue == "REFERENCE" && keepRef) || (typeValue == "NON-REFERENCE" && keepNonref)){
                            match2 = true;
                        }
                    }
                } else {
                    if (txtValue.toUpperCase().indexOf(filter1) > -1 || filter1 == ""){
                        if ((typeValue == "REFERENCE" && keepRef) || (typeValue == "NON-REFERENCE" && keepNonref)){
                            match1 = true;
                        }
                    }
                    if (txtValue.toUpperCase().indexOf(filter2) > -1 || filter2 == ""){
                        if ((typeValue == "REFERENCE" && keepRef) || (typeValue == "NON-REFERENCE" && keepNonref)){
                            match2 = true;
                        }
                    }
                }
            }
        }
        // console.log(match1);
        // console.log("filtered");
        if ((filterOption == "or" && (match1 || match2)) || (filterOption == "and" && (match1 && match2))){
            valuesArray[rows[1].getElementsByTagName("td").length] = "";
            rowArray[i] = valuesArray;

        } else {
            valuesArray[rows[1].getElementsByTagName("td").length] = "none";
            rowArray[i] = valuesArray;
        }
    }

    var t1 = performance.now();
    console.log("change display: " + (t1-t0) + "ms");

    // clear table
    for (i = rows.length-1; i > 0; i--){
        table.deleteRow(i);
    }

    var t2 = performance.now();
    console.log("clear table: " + (t2-t1) + "ms");

    // repolulate table with sorted values
    for (i = 1; i < rowArray.length; i++){
        var newRow = document.createElement("tr");
        var newCell;

        for (j = 0; j < rowArray[1].length-1; j++){
            newCell = newRow.insertCell(j);
            newCell.innerHTML = rowArray[i][j];
            newCell.classList.add("values");
        }
        newRow.style.display = rowArray[i][rowArray[1].length-1];
        table.tBodies[0].appendChild(newRow);
    }

    var t3 = performance.now();
    console.log("repopulate display: " + (t3-t2) + "ms");
    
    table.style.display = "";

    console.log("total duration: " + (t3-t0) + "ms");
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