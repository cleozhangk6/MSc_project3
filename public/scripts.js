
function populate_uniprot(uni){
    document.getElementById('input-uniprot').value=uni;
}

function toggleInfo() {
    var x = document.getElementById("info-box");
    if (x.style.visibility === "" || x.style.visibility === "hidden") {
        x.style.visibility = "visible";
    } else {
        x.style.visibility = "hidden";
    }
}

// Change color code accordingly
var colcodes = document.getElementsByClassName("colcode")
for (var i = 0; i < colcodes.length; i++) {
    colcodes.item(i).style.color = colcodes.item(i).textContent
}