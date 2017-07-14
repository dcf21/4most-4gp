// graphs_custom.js
// $Id$

// Library to display custom ephemeris graphs

function GraphCustom(source, context) {
    var self = this;

    this._source = source;
    this._context = context;
    this._width = 0;
    this._graph_data = [];

    $.getJSON(this._source, function (data) {
        self._graph_data = [];
        $.each(data, function (i,pt) {
            self._graph_data.push({'x': pt[0], 'y': pt[1]})
        });
        self._resize();
    });

    $(window).resize(function () {
        self._resize();
    });
}

GraphCustom.prototype._redraw = function () {
    var ctx = $("canvas", this._context);
    this._myLineChart = new Chart(ctx).Scatter({
        type: 'line',
        data: [{'data': this._graph_data}],
        options: {
            scales: {
                xAxes: [{
                    type: 'linear',
                    position: 'bottom',
                    scaleLabel: {
                        'display': true,
                        'labelString': "Wavelength / A"
                    }
                }],
                yAxes: [{
                    type: 'linear',
                    scaleLabel: {
                        'display': true,
                        'labelString': "Value"
                    }
                }]
            }
        }
    });
};

GraphCustom.prototype._resize = function () {
    this._width = this._context.width() - 20;
    $("canvas", this._context).width(this._width).height(300);
    this._redraw();
};


$(function () {
    $(".show_spectrum").each(function (i, el) {
        var elj = $(el);
        var meta = elj.data("source");
        var handler = new GraphCustom(meta, elj);
        elj.data("handler", handler);
    });
});

