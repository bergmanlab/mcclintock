function fillTable(data, tableId, maxTableSize){
    //  removes undefined rows
    var rawData = data.filter(function(x){ return x !== undefined;})

    console.log("filltable length:",rawData.length)
    console.log(rawData)


    var table = document.getElementById(tableId);
    maxRows = maxTableSize;
    if (maxTableSize > rawData.length-1){
        maxRows = rawData.length-1;
    }
    for (var i = 1; i <= maxRows; i++){
        var newRow = document.createElement("tr");
        for (let j = 0; j < rawData[0].length; j++){
            // var text = document.createTextNode(rawData[i][j]);
            var newCol = document.createElement('td');
            newCol.innerHTML = rawData[i][j];
            newCol.classList.add("values");
            // newCol.appendChild(text);
            newRow.appendChild(newCol);
        }
        table.appendChild(newRow);
    }
}

function sortTableLarge(rawData, tableId, colToSort, maxTableSize, btnPrefix, numeric=false){
    console.log(tableId,colToSort)
    var t0 = performance.now();
    var table = document.getElementById(tableId);
    var rows = table.getElementsByTagName("tr");

    // remove sorted label from non-sorted column headers
    for (var c = 0; c < rows[0].length; c++){
        if (c != colToSort){
            rows[0].getElementsByTagName("th")[c].classList.remove("sorted");
            rows[0].getElementsByTagName("th")[c].classList.remove("invSorted");
        }
    }

    //clear table
    for (i = rows.length-1; i > 0; i--){
        table.deleteRow(i);
    }

    // remove header row from sorting
    var rowArray = [...rawData];
    rowArray.shift();

    // sort data
    if (rows[0].getElementsByTagName("th")[colToSort].classList.contains("sorted")){
        console.log("is already sorted");
        sortByColumn(rowArray, colToSort, numeric=numeric, invert=true);
        rows[0].getElementsByTagName("th")[colToSort].classList.remove("sorted");
        rows[0].getElementsByTagName("th")[colToSort].classList.add("invSorted");
    } else {
        sortByColumn(rowArray, colToSort, numeric=numeric, invert=false);
        rows[0].getElementsByTagName("th")[colToSort].classList.add("sorted");
    }

    // fill table
    rowArray.unshift(rawData[0]);
    fillTable(rowArray, tableId, maxTableSize);

    // add button coloring to shade current section
    setupSectionButtons(rowArray, maxTableSize, btnPrefix, section=1);

    return rowArray;
}

function filterData(rawData, tableId, maxTableSize, btnPrefix, filterId1, filterId2, exactId, refId, nonrefId, singleField, containsBothId, predTypeRow, keepAllId){
    var filteredData = [];
    console.log(filterId1);
    var filter1 = document.getElementById(filterId1).value.toUpperCase();
    console.log(document.getElementById(filterId1).classList)
    var filter2 = document.getElementById(filterId2).value.toUpperCase();
    var exactMatch = document.getElementById(exactId).checked;
    var keepRef = document.getElementById(refId).checked;
    var keepNonref = document.getElementById(nonrefId).checked;

    if (containsBothId == false){
        containsBoth = false;
    } else {
        var containsBothValue = document.getElementById(containsBothId).value;
        var containsBoth = false;
        if (containsBothValue == "and"){
            containsBoth = true;
        }
    }

    if (keepAllId == ""){
       var keepAll = false;
    } else {
        var keepAll = document.getElementById(keepAllId).checked;
    }

    if (singleField){
        for (var i = 1; i < rawData.length; i++){
            var match =false;
            for (var j = 0; j < 2; j++){
                var valToCheck = rawData[i][j].toUpperCase();
                var predType = rawData[i][predTypeRow].toUpperCase();
                if (exactMatch){
                    if ((valToCheck == filter1) || (filter1 == "")){
                        if (((predType == "REFERENCE") && keepRef) || ((predType == "NON-REFERENCE") && keepNonref) || ((predType == "ALL") && keepAll)){
                            match = true;
                        }
                    }

                } else {
                    if ((valToCheck.indexOf(filter1) > -1) || (filter1 == "")){
                        if (((predType == "REFERENCE") && keepRef) || ((predType == "NON-REFERENCE") && keepNonref) || ((predType == "ALL") && keepAll)){
                            console.log("here:", filter1, valToCheck);
                            match = true;
                        }
                    }
                }
            }
            if (match){
                filteredData.push(rawData[i]);
            }
        }
    } else {
        for (var i = 1; i < rawData.length; i++){
            var match1 = false;
            var match2 = false;
            for (var j = 0; j < 2; j++){
                var valToCheck = rawData[i][j].toUpperCase();
                var predType = rawData[i][predTypeRow].toUpperCase();
                console.log(valToCheck, filter1, filter2);
                if (exactMatch){
                    if ((valToCheck == filter1) || (filter1 == "")){
                        if (((predType == "REFERENCE") && keepRef) || ((predType == "NON-REFERENCE") && keepNonref) || ((predType == "ALL") && keepAll)){
                            match1 = true;
                        }
                    }
                    if ((valToCheck == filter2) || (filter2 == "")){
                        if (((predType == "REFERENCE") && keepRef) || ((predType == "NON-REFERENCE") && keepNonref) || ((predType == "ALL") && keepAll)){
                            match2 = true;
                        }
                    }
                } else {
                    if ((valToCheck.indexOf(filter1) > -1) || (filter1 == "")){
                        if (((predType == "REFERENCE") && keepRef) || ((predType == "NON-REFERENCE") && keepNonref) || ((predType == "ALL") && keepAll)){
                            match1 = true;
                        }
                    }
                    if ((valToCheck.indexOf(filter2) > -1) || (filter2 == "")){
                        if (((predType == "REFERENCE") && keepRef) || ((predType == "NON-REFERENCE") && keepNonref) || ((predType == "ALL") && keepAll)){
                            match2 = true;
                        }
                    }
                }
            }
            console.log(match1, match2);
            if (containsBoth){
                if (match1 && match2){
                    filteredData.push(rawData[i]);
                }
            } else {
                if (match1 || match2){
                    filteredData.push(rawData[i]);
                }
            }
        }
    }

    filteredData.unshift(rawData[0]);
    console.log(filteredData.length);

    var table = document.getElementById(tableId);
    var rows = table.getElementsByTagName("tr");
    //clear table
    for (i = rows.length-1; i > 0; i--){
        table.deleteRow(i);
    }
    fillTable(filteredData, tableId, maxTableSize);
    console.log(filteredData.length)

    // add button coloring to shade current section
    setupSectionButtons(filteredData, maxTableSize, btnPrefix, section=1);

    return filteredData;
}


function showSection(inData, tableId, maxTableSize, btnId, btnPrefix){
    var dataToShow = [];
    var button = document.getElementById(btnId);
    var section = parseInt(button.value);

    var data = [...inData];
    data.shift();

    var metrics = calcSectionMetrics(data, maxTableSize);
    var start = metrics[0];
    var end = metrics[1];
    var binNum = metrics[2];

    for (var b = 1; b <= binNum; b++){
        if (b == section){
            for (x = start; x < end; x++){
                dataToShow.push(data[x]);
            }
        }
        start += maxTableSize;
        end += maxTableSize;
    }

    dataToShow.unshift(data[0]);

    var table = document.getElementById(tableId);
    var rows = table.getElementsByTagName("tr");
    //clear table
    for (i = rows.length-1; i > 0; i--){
        table.deleteRow(i);
    }

    fillTable(dataToShow, tableId, maxTableSize);

    // add button coloring to shade current section
    setupSectionButtons(inData, maxTableSize, btnPrefix, section=section);
}

function calcSectionMetrics(data, maxTableSize){
    var binNum = parseInt(data.length / maxTableSize);
    if ((data.length / maxTableSize) - binNum > 0){
        binNum += 1;
    }

    var start = 0;
    var end = maxTableSize;
    if (maxTableSize > data.length){
        end = data.length;
    }
    console.log(start, end, binNum);
    return [start, end, binNum];
}

function setupSectionButtons(inData, maxTableSize, btnPrefix, section=1){
    
    var data = [...inData];
    data.shift();
    var metrics = calcSectionMetrics(data, maxTableSize);
    var binNum = metrics[2];
    var hide1 = document.getElementById(btnPrefix+"h1");
    var hide2 = document.getElementById(btnPrefix+"h2");

    console.log("section:", section, "bins:", binNum);

    // remove current page styling
    for (var x = 1; x < 8; x++){
        var button = document.getElementById(btnPrefix+x.toString());
        button.classList.remove('currentPage');
    }

    // no hidden sections
    if (binNum < 8){
        hide1.style.display = "none";
        hide2.style.display = "none";
        for (var x = 1; x < 8; x++){
            var button = document.getElementById(btnPrefix+x.toString());
            if (x == section){
                button.classList.add('currentPage');
            }
            if (x > binNum){
                button.style.display = "none";
            }
        }
    } else {
        hide1.style.display = "";
        hide2.style.display = "";

        // hide downstream sections
        if (section < 5){
            hide1.style.display = "none";
            for (var x = 1; x < 7; x++){
                setButton(btnPrefix+x.toString(), x);
                if (x == section){
                    document.getElementById(btnPrefix+x.toString()).classList.add("currentPage");
                }
            }
            setButton(btnPrefix+"7", binNum);

        // hide upstream sections
        } else if (section > (binNum-4)) {
            hide2.style.display = "none";
            setButton(btnPrefix+"1", 1);

            for (var x = 2; x < 8; x++){
                setButton(btnPrefix+x.toString(), (binNum-(7-x)));
                if ((binNum-(7-x)) == section){
                    document.getElementById(btnPrefix+x.toString()).classList.add("currentPage");
                }
            }
        
        // hide up and downstream sections
        } else {
            document.getElementById(btnPrefix+"4").classList.add("currentPage");
            setButton(btnPrefix+"1", 1);
            setButton(btnPrefix+"2", section-2);
            setButton(btnPrefix+"3", section-1);
            setButton(btnPrefix+"4", section);
            setButton(btnPrefix+"5", section+1);
            setButton(btnPrefix+"6", section+2);
            setButton(btnPrefix+"7", binNum);
        }
    }
}

function setButton(buttonId, value){
    document.getElementById(buttonId).innerHTML = value;
    document.getElementById(buttonId).value = value;
}


function isSorted(arr, col, numeric=false, inverted=false){
    for (var i = 0; i < arr.length -1; i++){
        if (numeric){
            if (inverted){
                if (parseInt(arr[i][col]) < parseInt(arr[i+1][col])){
                    console.log(arr[i][col], arr[i+1][col])
                    return false;
                }
            } else {
                if (parseInt(arr[i][col]) > parseInt(arr[i+1][col])){
                    console.log(arr[i][col], arr[i+1][col])
                    return false;
                }
            }

        } else {
            if (inverted){
                if (arr[i][col].toLowerCase() < arr[i+1][col].toLowerCase()){
                    console.log(arr[i][col], arr[i+1][col])
                    return false;
                }
            } else {
                if (arr[i][col].toLowerCase() > arr[i+1][col].toLowerCase()){
                    console.log(arr[i][col], arr[i+1][col])
                    return false;
                }
            }

        }
    }
    return true;
}


function sortTable(id, n, numeric=false){
    console.log(id,n)
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

        table.appendChild(newRow);
    }

    var t5 = performance.now();
    console.log("write table: " + (t5-t4) + "ms");

    table.style.display = "";
    var t7 = performance.now();
    console.log("total duration: "+ (t7-t0) + "ms");
}


function invertArray(arr){
    var invertedArray = [...arr];
    for (var i = 0; i < arr.length; i++){
        invertedArray[i] = arr[(arr.length)-i]
    }

    return invertedArray;
}

function sortByColumn(arr, colIdx, numeric=false, invert=false){
    if (invert){
        if (numeric){
            function sortFunction(a,b){
                if (parseInt(a[colIdx]) < parseInt(b[colIdx])){
                    return 1;
                } else {
                    return -1;
                }
            }
        } else {
            function sortFunction(a,b){
                if (a[colIdx].toLowerCase() < b[colIdx].toLowerCase()){
                    return 1;
                } else {
                    return -1;
                }
            }
        }
    } else {
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
    }

    arr.sort(sortFunction);
    // tmp = arr[-1];
    // arr.unshift(tmp);
    // arr.pop();
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
        table.appendChild(newRow);
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