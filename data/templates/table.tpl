{# Update the template_structure.html document too #}
{%- block before_style -%}{%- endblock before_style -%}
{% block style %}
<link rel="stylesheet" href="https://fonts.googleapis.com/css?family={{font.replace(" ", "+")}}"/>
<style  type="text/css" >
/* line 2, ../sass/_sortable.sass */
table[data-sortable] {
  font-family: '{{font}}';
  font-size: 90%;
  border-collapse: collapse;
  border-spacing: 0;
}
/* line 6, ../sass/_sortable.sass */
table[data-sortable] th {
  vertical-align: bottom;
  font-weight: bold;
}
/* line 10, ../sass/_sortable.sass */
table[data-sortable] th, table[data-sortable] td {
  text-align: left;
  padding: 2px;
}
/* line 14, ../sass/_sortable.sass */
table[data-sortable] th:not([data-sortable="false"]) {
  -webkit-user-select: none;
  -moz-user-select: none;
  -ms-user-select: none;
  -o-user-select: none;
  user-select: none;
  -webkit-tap-highlight-color: rgba(0, 0, 0, 0);
  -webkit-touch-callout: none;
  cursor: pointer;
}
/* line 26, ../sass/_sortable.sass */
table[data-sortable] th:after {
  content: "";
  visibility: hidden;
  display: inline-block;
  vertical-align: inherit;
  height: 0;
  width: 0;
  border-width: 5px;
  border-style: solid;
  border-color: transparent;
  margin-right: 1px;
  margin-left: 10px;
  float: right;
}
/* line 40, ../sass/_sortable.sass */
table[data-sortable] th[data-sorted="true"]:after {
  visibility: visible;
}
/* line 43, ../sass/_sortable.sass */
table[data-sortable] th[data-sorted-direction="descending"]:after {
  border-top-color: inherit;
  margin-top: 8px;
}
/* line 47, ../sass/_sortable.sass */
table[data-sortable] th[data-sorted-direction="ascending"]:after {
  border-bottom-color: inherit;
  margin-top: 3px;
}

/* line 6, ../sass/sortable-theme-slick.sass */
table[data-sortable].sortable-theme-slick {
  color: #333333;
  background: white;
  border: 1px solid #e0e0e0;
}
/* line 11, ../sass/sortable-theme-slick.sass */
table[data-sortable].sortable-theme-slick thead th {
  background-image: -webkit-gradient(linear, 50% 0%, 50% 100%, color-stop(0%, #ffffff), color-stop(100%, #eeeeee));
  background-image: -webkit-linear-gradient(#ffffff, #eeeeee);
  background-image: -moz-linear-gradient(#ffffff, #eeeeee);
  background-image: -o-linear-gradient(#ffffff, #eeeeee);
  background-image: linear-gradient(#ffffff, #eeeeee);
  background-color: #f0f0f0;
}
/* line 16, ../sass/sortable-theme-slick.sass */
table[data-sortable].sortable-theme-slick tbody td {
  border-top: 1px solid #e0e0e0;
}
/* line 19, ../sass/sortable-theme-slick.sass */
table[data-sortable].sortable-theme-slick tbody > tr:nth-child(odd) > td {
  background-color: #f9f9f9;
}
/* line 22, ../sass/sortable-theme-slick.sass */
table[data-sortable].sortable-theme-slick th[data-sorted="true"] {
  -webkit-box-shadow: inset 1px 0 #bce8f1, inset -1px 0 #bce8f1;
  -moz-box-shadow: inset 1px 0 #bce8f1, inset -1px 0 #bce8f1;
  box-shadow: inset 1px 0 #bce8f1, inset -1px 0 #bce8f1;
  color: #3a87ad;
  background: #d9edf7;
  border-bottom-color: #bce8f1;
}
/* line 28, ../sass/sortable-theme-slick.sass */
table[data-sortable].sortable-theme-slick th[data-sorted="true"]:first-child {
  -webkit-box-shadow: inset -1px 0 #bce8f1;
  -moz-box-shadow: inset -1px 0 #bce8f1;
  box-shadow: inset -1px 0 #bce8f1;
}
/* line 31, ../sass/sortable-theme-slick.sass */
table[data-sortable].sortable-theme-slick th[data-sorted="true"]:last-child {
  -webkit-box-shadow: inset 1px 0 #bce8f1;
  -moz-box-shadow: inset 1px 0 #bce8f1;
  box-shadow: inset 1px 0 #bce8f1;
}
/* line 34, ../sass/sortable-theme-slick.sass */
table[data-sortable].sortable-theme-slick th[data-sorted="true"][data-sorted-direction="descending"]:after {
  border-top-color: #3a87ad;
}
/* line 37, ../sass/sortable-theme-slick.sass */
table[data-sortable].sortable-theme-slick th[data-sorted="true"][data-sorted-direction="ascending"]:after {
  border-bottom-color: #3a87ad;
}


{% block col_heading_style %}
    .{{col_heading_style.name}} {
    {% for p,val in col_heading_style.props %}
      {{p}}: {{val}};
    {% endfor -%}
    }
{% endblock col_heading_style %}

{% block circle_styles %}
{% for s in circle_styles %}
    .{{s.name}} {
    {% for p,val in s.props %}
      {{p}}: {{val}};
    {% endfor -%}
    }
{%- endfor -%}
{% endblock circle_styles %}
{% block palette_styles %}
{% for s in palette_styles %}
    .{{s.name}} {
    {% for p,val in s.props %}
      {{p}}: {{val}};
    {% endfor -%}
    }
{%- endfor -%}
{% endblock palette_styles %}


{% block table_styles %}
{% for s in table_styles %}
    #T_{{uuid}} {{s.selector}} {
    {% for p,val in s.props %}
      {{p}}: {{val}};
    {% endfor -%}
    }
{%- endfor -%}
{% endblock table_styles %}
{% block before_cellstyle %}{% endblock before_cellstyle %}
{% block cellstyle %}
{%- for s in cellstyle %}
    {%- for selector in s.selectors -%}{%- if not loop.first -%},{%- endif -%}#T_{{uuid}}{{selector}}{%- endfor -%} {
    {% for p,val in s.props %}
        {{p}}: {{val}};
    {% endfor %}
    }
{%- endfor -%}
{%- endblock cellstyle %}
</style>
{%- endblock style %}
{%- block before_table %}{% endblock before_table %}
{%- block table %}
<table id="T_{{uuid}}" {% if table_attributes %}{{ table_attributes }}{% endif %}>
{%- block caption %}
{%- if caption -%}
    <caption>{{caption}}</caption>
{%- endif -%}
{%- endblock caption %}
{%- block thead %}
<thead>
    {%- block before_head_rows %}{% endblock %}
    {%- for r in head %}
    {%- block head_tr scoped %}
    <tr>
        {%- for c in r %}
        {%- if c.is_visible != False %}
        <{{ c.type }} class="{{c.class}}" {{ c.attributes|join(" ") }}>{{c.value}}</{{ c.type }}>
        {%- endif %}
        {%- endfor %}
    </tr>
    {%- endblock head_tr %}
    {%- endfor %}
    {%- block after_head_rows %}{% endblock %}
</thead>
{%- endblock thead %}
{%- block tbody %}
<tbody>
    {% block before_rows %}{% endblock before_rows %}
    {% for r in body %}
    {% block tr scoped %}
    <tr>
        {% for c in r %}
        {% if c.is_visible != False %}
        <{{ c.type }} {% if c.id is defined -%} id="T_{{ uuid }}{{ c.id }}" {%- endif %} class="{{ c.class }}" {{ c.attributes|join(" ") }}>{{ c.display_value }}</{{ c.type }}>
        {% endif %}
        {%- endfor %}
    </tr>
    {% endblock tr %}
    {%- endfor %}
    {%- block after_rows %}{%- endblock after_rows %}
</tbody>
{%- endblock tbody %}
</table>
{%- endblock table %}
{%- block after_table %}{% endblock after_table %}
<script>
/*! sortable.js 0.8.0 */
(function(){var a,b,c,d,e,f,g;a="table[data-sortable]",d=/^(-?[£$¤]?[\d,.e\-]+%?|inf)$/,g=/^\s+|\s+$/g,c=["click"],f="ontouchstart"in document.documentElement,f&&c.push("touchstart"),b=function(a,b,c){return null!=a.addEventListener?a.addEventListener(b,c,!1):a.attachEvent("on"+b,c)},e={init:function(b){var c,d,f,g,h;for(null==b&&(b={}),null==b.selector&&(b.selector=a),d=document.querySelectorAll(b.selector),h=[],f=0,g=d.length;g>f;f++)c=d[f],h.push(e.initTable(c));return h},initTable:function(a){var b,c,d,f,g,h;if(1===(null!=(h=a.tHead)?h.rows.length:void 0)&&"true"!==a.getAttribute("data-sortable-initialized")){for(a.setAttribute("data-sortable-initialized","true"),d=a.querySelectorAll("th"),b=f=0,g=d.length;g>f;b=++f)c=d[b],"false"!==c.getAttribute("data-sortable")&&e.setupClickableTH(a,c,b);return a}},setupClickableTH:function(a,d,f){var g,h,i,j,k,l;for(i=e.getColumnType(a,f),h=function(b){var c,g,h,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z,A,B,C,D;if(b.handled===!0)return!1;for(b.handled=!0,m="true"===this.getAttribute("data-sorted"),n=this.getAttribute("data-sorted-direction"),h=m?"ascending"===n?"descending":"ascending":i.defaultSortDirection,p=this.parentNode.querySelectorAll("th"),s=0,w=p.length;w>s;s++)d=p[s],d.setAttribute("data-sorted","false"),d.removeAttribute("data-sorted-direction");if(this.setAttribute("data-sorted","true"),this.setAttribute("data-sorted-direction",h),o=a.tBodies[0],l=[],m){for(D=o.rows,v=0,z=D.length;z>v;v++)g=D[v],l.push(g);for(l.reverse(),B=0,A=l.length;A>B;B++)k=l[B],o.appendChild(k)}else{for(r=null!=i.compare?i.compare:function(a,b){return b-a},c=function(a,b){return a[0]===b[0]?a[2]-b[2]:i.reverse?r(b[0],a[0]):r(a[0],b[0])},C=o.rows,j=t=0,x=C.length;x>t;j=++t)k=C[j],q=e.getNodeValue(k.cells[f]),null!=i.comparator&&(q=i.comparator(q)),l.push([q,k,j]);for(l.sort(c),u=0,y=l.length;y>u;u++)k=l[u],o.appendChild(k[1])}return"function"==typeof window.CustomEvent&&"function"==typeof a.dispatchEvent?a.dispatchEvent(new CustomEvent("Sortable.sorted",{bubbles:!0})):void 0},l=[],j=0,k=c.length;k>j;j++)g=c[j],l.push(b(d,g,h));return l},getColumnType:function(a,b){var c,d,f,g,h,i,j,k,l,m,n;if(d=null!=(l=a.querySelectorAll("th")[b])?l.getAttribute("data-sortable-type"):void 0,null!=d)return e.typesObject[d];for(m=a.tBodies[0].rows,h=0,j=m.length;j>h;h++)for(c=m[h],f=e.getNodeValue(c.cells[b]),n=e.types,i=0,k=n.length;k>i;i++)if(g=n[i],g.match(f))return g;return e.typesObject.alpha},getNodeValue:function(a){var b;return a?(b=a.getAttribute("data-value"),null!==b?b:"undefined"!=typeof a.innerText?a.innerText.replace(g,""):a.textContent.replace(g,"")):""},setupTypes:function(a){var b,c,d,f;for(e.types=a,e.typesObject={},f=[],c=0,d=a.length;d>c;c++)b=a[c],f.push(e.typesObject[b.name]=b);return f}},e.setupTypes([{name:"numeric",defaultSortDirection:"descending",match:function(a){return a.match(d)},comparator:function(a){if(a=="inf"){return Infinity}else{return parseFloat(a.replace(/^[^0-9\-]+/g,""),10)||0}}},{name:"date",defaultSortDirection:"ascending",reverse:!0,match:function(a){return!isNaN(Date.parse(a))},comparator:function(a){return Date.parse(a)||0}},{name:"alpha",defaultSortDirection:"ascending",match:function(){return!0},compare:function(a,b){return a.localeCompare(b)}}]),setTimeout(e.init,0),"function"==typeof define&&define.amd?define(function(){return e}):"undefined"!=typeof exports?module.exports=e:window.Sortable=e}).call(this);
</script>
