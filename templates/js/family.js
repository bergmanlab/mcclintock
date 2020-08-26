function sortTable(id, n, numeric=false){
    console.log(id)
    var table = document.getElementById(id);
    var swapped = true;
    var sortDirection = "lowToHigh";

    var rows,i, shouldSwap,x,y, switchCount = 0;
    while(swapped){
        swapped = false;
        rows = table.rows;
        for (i =1; i < (rows.length-1); i++){
            shouldSwap = false;

            x = rows[i].getElementsByTagName("TD")[n];
            y = rows[i + 1].getElementsByTagName("TD")[n];

            if (sortDirection == "lowToHigh"){
                if (numeric){
                    if (parseInt(x.innerHTML) > parseInt(y.innerHTML)){
                        shouldSwap = true;
                        break;
                    }
                } else {
                    if (x.innerHTML.toLowerCase() > y.innerHTML.toLowerCase()){
                        shouldSwap = true;
                        break;
                    }
                }

            } else if (sortDirection == "highToLow"){
                if (numeric){
                    if (parseInt(x.innerHTML) < parseInt(y.innerHTML)){
                        shouldSwap = true;
                        break;
                    }
                } else {
                    if (x.innerHTML.toLowerCase() < y.innerHTML.toLowerCase()){
                        shouldSwap = true;
                        break;
                    }
                }
            }
        }
        if (shouldSwap){
            rows[i].parentNode.insertBefore(rows[i + 1], rows[i]);
            swapped = true;
            switchCount++;
        } else {
            console.log(switchCount)
            console.log(sortDirection)
            if (switchCount == 0 && sortDirection == "lowToHigh"){
                sortDirection = "highToLow";
                swapped = true;
                console.log(sortDirection);
            }
        }
    }
}

function filterTable(inputId, tableId, checkbox){
    var input = document.getElementById(inputId);
    var filter = input.value.toUpperCase();
    var table = document.getElementById(tableId);
    var rows = table.getElementsByTagName("tr");
    var exactMatch = document.getElementById(checkbox).checked;
    console.log(exactMatch);
    var display = false;
    var data;
    for (var i = 1; i < rows.length; i++){
        display = false;
        for (var j = 0; j < rows[i].getElementsByTagName("td").length; j++){
            data = rows[i].getElementsByTagName("td")[j];
            if (data){
                txtValue = data.textContent || data.innerText;
                if (exactMatch){
                    if (txtValue.toUpperCase() == filter){
                        display = true;
                    }
                } else {
                    if (txtValue.toUpperCase().indexOf(filter) > -1){
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