function fillTable(data, tableId, maxTableSize){
    //  removes undefined rows
    var rawData = data.filter(function(x){ return x !== undefined;})
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
    var filter1 = document.getElementById(filterId1).value.toUpperCase();
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

    var table = document.getElementById(tableId);
    var rows = table.getElementsByTagName("tr");
    //clear table
    for (i = rows.length-1; i > 0; i--){
        table.deleteRow(i);
    }
    fillTable(filteredData, tableId, maxTableSize);

    // add button coloring to shade current section
    setupSectionButtons(filteredData, maxTableSize, btnPrefix, section=1);

    return filteredData;
}


function showSection(inData, tableId, maxTableSize, btnId, btnPrefix){
    var dataToShow = [];
    if (!isNaN(btnId)){
        var section = parseInt(btnId);
    } else {
        var button = document.getElementById(btnId);
        var section = parseInt(button.value);
    }


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
    console.log(inData.length, maxTableSize, btnPrefix, section);
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

    return [start, end, binNum];
}

function setupSectionButtons(inData, maxTableSize, btnPrefix, section=1){
    
    var data = [...inData];
    data.shift();
    var metrics = calcSectionMetrics(data, maxTableSize);
    var binNum = metrics[2];
    var hide1 = document.getElementById(btnPrefix+"h1");
    var hide2 = document.getElementById(btnPrefix+"h2");

    // remove current page styling; reset hidden buttons
    for (var x = 1; x < 8; x++){
        var button = document.getElementById(btnPrefix+x.toString());
        button.classList.remove('currentPage');
        button.style.display = "";
    }

    // resets hidden input and input submit
    var button = document.getElementById(btnPrefix+"Go");
    button.style.display = "";
    var button = document.getElementById(btnPrefix+"Input");
    button.style.display = "";

    // no hidden sections
    if (binNum < 8){
        hide1.style.display = "none";
        hide2.style.display = "none";
        for (var x = 1; x < 8; x++){
            var button = document.getElementById(btnPrefix+x.toString());
            setButton(btnPrefix+x.toString(), x);
            button.style.display = "";
            if (x == section){
                button.classList.add('currentPage');
            }
            if (x > binNum){
                button.style.display = "none";
            }
            // if only one section, hide the section button, submit button, and input
            if (binNum == 1 && x == 1){
                button.style.display = "none";
                document.getElementById(btnPrefix+"Go").style.display = "none";
                document.getElementById(btnPrefix+"Input").style.display = "none";
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

function goToSection(inputId, inData, table, maxTableSize, btnPrefix){
    console.log("goToSection");
    var input  = document.getElementById(inputId);
    var section = input.value;
    var data = [...inData];
    data.shift();

    input.classList.remove("error");
    var metrics = calcSectionMetrics(data, maxTableSize);
    var start = metrics[0];
    var end = metrics[1];
    var binNum = metrics[2];

    if (!isNaN(section) && section <= binNum && section > 0){
        showSection(inData, table, maxTableSize, section, btnPrefix);
    } else {
        console.log("Error:", section);
        input.classList.add("error");
    }
}