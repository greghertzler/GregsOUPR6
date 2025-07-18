// Javascript for Topics

window.addEventListener('message', TopicReceiveFindHelp, false);

function TopicReceiveFindHelp(event) {
    if (typeof event.data === 'string') {
        var who = event.data.slice(0, 3);
        var what = event.data.slice(3);
        if (who == 'fnd') { this.TopicFind(what); }
        if (who == 'hlp') {
            if (what == 'btnclk') { this.SendFocus(); }
            else { this.TopicOpenTo(what); }
        }
    }
}

function SendTopic() {
    var msg = 'tpc' + getTopicId();
    window.parent.frames[0].postMessage(msg, '*');
    window.parent.frames[1].postMessage(msg, '*');
}

function SendFocus() {
    var msg = 'tpc' + getTopicId();
    window.parent.frames[1].postMessage(msg, '*');
}

function TopicOpenTo(tpcid) {
    window.location.replace(tpcid + '.html');
}

function TopicFind(FindText) {
    var tpcp = document.getElementsByTagName('p');
    var tpcli = document.getElementsByTagName('li');
    var tpccap = document.getElementsByTagName('caption');
    var tpcth = document.getElementsByTagName('th');
    var tpctd = document.getElementsByTagName('td');
    var tpcdt = document.getElementsByTagName('dt');
    var tpcdd = document.getElementsByTagName('dd');
    var tpcasi = document.getElementsByTagName('aside');
    var i;
    for (i=0; i<tpcp.length; i++) {
        tpcp[i].setAttribute('style', 'border-style: none');
    }
    for (i=0; i<tpcli.length; i++) {
        tpcli[i].setAttribute('style', 'border-style: none');
    }
    for (i=0; i<tpccap.length; i++) {
        tpccap[i].setAttribute('style', 'border-style: none');
    }
    for (i=0; i<tpcth.length; i++) {
        tpcth[i].setAttribute('style', 'border-style: none');
    }
    for (i=0; i<tpctd.length; i++) {
        tpctd[i].setAttribute('style', 'border-style: none');
    }
    for (i = 0; i < tpcdt.length; i++) {
        tpcdt[i].setAttribute('style', 'border-style: none');
    }
    for (i = 0; i < tpcdd.length; i++) {
        tpcdd[i].setAttribute('style', 'border-style: none');
    }
    for (i = 0; i < tpcasi.length; i++) {
        tpcasi[i].setAttribute('style', 'border-style: none');
    }
    if (FindText != '') {
        var str;
        var n;
        for (i=0; i<tpcp.length; i++) {
            str = tpcp[i].textContent.toLowerCase();
            n = str.search(FindText.toLowerCase());
            if (n>-1) {
                tpcp[i].setAttribute('style', 'border-style: solid');
            }
        }
        for (i=0; i<tpcli.length; i++) {
            str = tpcli[i].textContent.toLowerCase();
            n = str.search(FindText.toLowerCase());
            if (n>-1) {
                tpcli[i].setAttribute('style', 'border-style: solid');
            }
        }
        for (i=0; i<tpcth.length; i++) {
            str = tpcth[i].textContent.toLowerCase();
            n = str.search(FindText.toLowerCase());
            if (n>-1) {
                tpcth[i].setAttribute('style', 'border-style: solid');
            }
        }
        for (i=0; i<tpccap.length; i++) {
            str = tpccap[i].textContent.toLowerCase();
            n = str.search(FindText.toLowerCase());
            if (n>-1) {
                tpccap[i].setAttribute('style', 'border-style: solid');
            }
        }
        for (i=0; i<tpctd.length; i++) {
            str = tpctd[i].textContent.toLowerCase();
            n = str.search(FindText.toLowerCase());
            if (n>-1) {
                tpctd[i].setAttribute('style', 'border-style: solid');
            }
        }
        for (i = 0; i < tpcdt.length; i++) {
            str = tpcdt[i].textContent.toLowerCase();
            n = str.search(FindText.toLowerCase());
            if (n > -1) {
                tpcdt[i].setAttribute('style', 'border-style: solid');
            }
        }
        for (i = 0; i < tpcdd.length; i++) {
            str = tpcdd[i].textContent.toLowerCase();
            n = str.search(FindText.toLowerCase());
            if (n > -1) {
                tpcdd[i].setAttribute('style', 'border-style: solid');
            }
        }
        for (i = 0; i < tpcasi.length; i++) {
            str = tpcasi[i].textContent.toLowerCase();
            n = str.search(FindText.toLowerCase());
            if (n > -1) {
                tpcasi[i].setAttribute('style', 'border-style: solid');
            }
        }
    }
}
