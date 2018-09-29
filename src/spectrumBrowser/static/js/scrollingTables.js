// scrollingTables.js
// $Id$

// This code implements scrolling tables

function ScrollingTable(context) {
    var self = this;
    this._context = context;
    this.header = $(".scrolltable_thead", context);

    this.body = $(".scrolltable_tbody", context);
    this.body_cols = $("td", $($("tbody tr", this.body)[0])); // Columns of first row of table's tbody

    this._refresh();

    $(window).resize(function () {
        self._refresh();
    });

    this.body.scroll(function () {
        self.scroll_table();
    });
}

ScrollingTable.prototype._refresh = function () {
    var self = this;

    // Display thead
    $("td", self._context).css("min-width", "inherit").css("max-width", "inherit");

    // Read widths of table columns
    self.column_widths = [];
    self.body_cols.each(function (i, el) {
        var width = $(el).outerWidth() + 12;
        var min_width = $(el).data("min-width");
        var max_width = $(el).data("max-width");
        if (min_width && (width < min_width)) width = min_width;
        if (max_width && (width > max_width)) width = max_width;
        self.column_widths.push(width);
    });

    // Update widths of header table columns to match body
    var rows = $("tr", self._context).not(".no_setwidth");
    $.each(self.column_widths, function (i, item) {
        $("td:nth-child(" + (i + 1) + ")",
            rows).css("overflow", "hidden").css("min-width", item).css("max-width", item);
    });

    self.scroll_table();
};

ScrollingTable.prototype.scroll_table = function () {

    // If user scrolls horizontally, we need to scroll the duplicate headers too

    // First set up a holder to span the entire width of the visible table
    this.header.css("left", this.body.offset().left).css("width", this.body.width()).css("overflow", "hidden");

    // Then offset the duplicate headers within this holder
    $("table", this.header).css("position", "relative").css("left", -this.body.scrollLeft());
};

// Initialise all HTML elements with class scrolltable
function scrollingTableRegister() {
    $(".scrolltable").each(function (i, el) {
        var elj = $(el);
        var handler = new ScrollingTable(elj);
        elj.data("handler", handler);
    });
}

$(scrollingTableRegister);
