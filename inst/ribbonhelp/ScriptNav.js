// Javascript for OUP_PaneNav

window.addEventListener('message', NavReceiveFindTopic, false);
window.addEventListener('click', NavExpandCollapse, false);
var lnks = document.getElementsByTagName('a');
var navp = document.getElementsByTagName('p');
var toc = [];
var undo = [];
var find;
var m;
var n;

window.onload = function () {
    var i;
    find = false;
    m = lnks.length;
    n = navp.length;
    toc[0] = 0; toc[1] = 0; toc[2] = 0; toc[3] = 1; toc[4] = 2; toc[5] = 2; toc[6] = 2; toc[7] = 3; toc[8] = 2; toc[9] = 3; toc[10] = 2; toc[11] = 2; toc[12] = 2; toc[13] = 2; toc[14] = 3; toc[15] = 3; toc[16] = 2; toc[17] = 2; toc[18] = 3; toc[19] = 3; toc[20] = 3; toc[21] = 3; toc[22] = 3; toc[23] = 1; toc[24] = 2; toc[25] = 2; toc[26] = 3; toc[27] = 3; toc[28] = 1; toc[29] = 2; toc[30] = 2; toc[31] = 2; toc[32] = 2; toc[33] = 1; toc[34] = 2; toc[35] = 3; toc[36] = 3; toc[37] = 2; toc[38] = 3; toc[39] = 2; toc[40] = 3; toc[41] = 1; toc[42] = 2; toc[43] = 2; toc[44] = 3; toc[45] = 2; toc[46] = 1; toc[47] = 2; toc[48] = 2; toc[49] = 2; toc[50] = 2; toc[51] = 2; toc[52] = 2; toc[53] = 2; toc[54] = 3; toc[55] = 3; toc[56] = 3; toc[57] = 3; toc[58] = 3; toc[59] = 3; toc[60] = 3;
    for (i = 0; i < m; i++) {
        if (toc[i] > 1) {
            lnks[i].style.display = 'none';
            undo[i] = 'none';
        }
        else {
            lnks[i].style.display = 'block';
            undo[i] = 'block';
        }
    }
}

function NavExpandCollapse(event) {
    if (!find) {
        var lnk = document.activeElement;
        var i;
        var hit = false;
        for (i = 0; hit == false & i < m; i++) {
            if (lnk.id == lnks[i].id) { hit = true; }
        }
        if (i < m) {
            var thisToc = toc[i - 1];
            var nextToc = toc[i];
            // hide or show
            if (thisToc < nextToc) {
                var j;
                hit = true;
                for (j = i; hit == true & toc[j] >= nextToc & j < m; j++) {
                    if (toc[j] == nextToc & lnks[j].style.display === 'none') { hit = false; }
                }
                // hide
                if (hit == true) {
                    for (j = i; toc[j] >= nextToc & j < m; j++) {
                        lnks[j].style.display = 'none';
                        undo[j] = 'none';
                    }
                }
                // show
                else {
                    for (j = i; toc[j] > thisToc & j < m; j++) {
                        if (toc[j] == nextToc) {
                            lnks[j].style.display = 'block';
                            undo[j] = 'block';
                        }
                    }
                }
            }
        }
    }
}

function NavReceiveFindTopic(event) {
    if (typeof event.data === 'string') {
        var who = event.data.slice(0, 3);
        var what = event.data.slice(3);
        if (who == 'fnd') { this.NavFind(what); }
        if (who == 'tpc') { this.NavFocus(what); }
    }
}

function NavFocus(tpcid) {
    var lnk = document.getElementById(tpcid);
    if (!find) {
        var i;
        var j;
        var hit = false;
        for (i = 0; hit == false & i < m; i++) {
            if (lnk.id == lnks[i].id) { hit = true; }
        }
        var thisToc = toc[i - 1];
        for (j = i - 1; toc[j] >= thisToc & j > 0; j--) { }
        for (i = j + 1; toc[i] >= thisToc & i < m; i++) {
            if (toc[i] == thisToc) {
                lnks[i].style.display = 'block';
                undo[i] = 'block';
            }
        }
    }
    lnk.focus();
    var lnkscr = window.pageYOffset;
    var lnktop = lnk.offsetTop;
    var lnkhgt = lnk.offsetHeight;
    var lnkprv = lnktop - lnkhgt;
    var lnknxt = lnktop + 2 * lnkhgt;
    var winscr = window.document.documentElement.scrollHeight;
    var maxscr = winscr - window.innerHeight;
    if (lnkprv > 0 && lnkprv < lnkscr) {
        window.scrollBy(0, -lnkhgt);
    }
    if (winscr - lnknxt > 0 && winscr - lnknxt < maxscr - lnkscr) {
        window.scrollBy(0, lnkhgt);
    }
}

function NavFind(findtext) {
    var i;
    var j;
    for (i = 0; i < n; i = i + 2) { navp[i].setAttribute('style', 'border-right-style: none'); }
    if (findtext != '') {
        var str;
        var k;
        var j = 0;
        for (i = 0; i < n; i = i + 2) {
            str = navp[i].textContent.toLowerCase();
            k = str.search(findtext.toLowerCase());
            if (k < 0) {
                str = navp[i + 1].textContent.toLowerCase();
                k = str.search(findtext.toLowerCase());
            }
            if (k > -1) {
                lnks[j].style.display = 'block';
                navp[i].setAttribute('style', 'border-right-style: solid');
            }
            else { lnks[j].style.display = 'none'; }
            j = j + 1;
        }
        find = true;
    }
    else {
        for (j = 0; j < m; j++) { lnks[j].style.display = undo[j]; }
        find = false;
    }
}
