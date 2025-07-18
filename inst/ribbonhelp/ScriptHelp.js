// Javascript for OUP_Help

var left;

window.onload = function ()
{
    left = 232;
    this.SendOpenTo();
    this.SetLeftWidth();
}

window.onresize = function ()
{
    this.SetLeftWidth();
}

window.onclick = function ()
{
    if (left == 232) { left = 0; }
    else {
        left = 232;
        this.SendFocus()
    }
    this.SetLeftWidth();
}

function ButtonClick()
{
    if (left == 232) { left = 0; }
    else {
        left = 232;
        this.SendFocus()
    }
    this.SetLeftWidth();
}

function SendFocus()
{
    var msg = 'hlpbtnclk';
    window.parent.frames[2].postMessage(msg, '*');
}

function SendOpenTo()
{
    var query = window.location.search;
    var tpcid = query.slice(1);
    if (tpcid != '') {
        var msg = 'hlp' + tpcid;
        window.parent.frames[2].postMessage(msg, '*');
    }
}

function SetLeftWidth()
{
    var width = window.innerWidth - left - 8;
    var tpcfrm = document.getElementById('tpc');
    var button = document.getElementById('btn');
    tpcfrm.style.left = left + 'px';
    tpcfrm.style.width = width + 'px';
    button.style.left = left + 'px';
    if (left == 0) { button.innerHTML = '&Gt;' }
    else { button.innerHTML = '&Lt;' }
}

