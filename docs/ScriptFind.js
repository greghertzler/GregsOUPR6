// Javascript for OUP_PaneFind

window.addEventListener('message', FindReceiveTopic, false);

function FindReceiveTopic(event) {
    if (typeof event.data === 'string') {
        var who = event.data.slice(0, 3);
        if (who == 'tpc') {
            var ft = document.getElementById('FindText');
            if (ft.value != '') { event.source.postMessage('fnd' + ft.value, '*'); }
        }
    }
}

function SendFind() {
    var ft = document.getElementById('FindText');
    var msg = 'fnd' + ft.value;
    window.parent.frames[1].postMessage(msg, '*');
    window.parent.frames[2].postMessage(msg, '*');
}

function FindFind() {

    this.SendFind();

    var ft = document.getElementById('FindText');
    if (ft.value == '') {
        ft.setAttribute('placeholder', 'Find...');
        ft.focus();
    }
}

function FindClear() {
    var ft = document.getElementById('FindText');
    ft.value = ''
    ft.setAttribute('placeholder', 'Find...');
    ft.focus();

    this.SendFind();

}
function KeyCheck(event) {
    var key = event.which;
    if (key == 13) { this.FindFind(); }
    if (key == 27) { this.FindClear(); }
}

function FindMouse(event) {
    if (event.offsetX > 150) {
        var ft = document.getElementById('FindText');
        if (ft.value != '') { ft.setAttribute('placeholder', '...Clear'); }
    }
}
