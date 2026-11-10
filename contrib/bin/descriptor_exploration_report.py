"""Offline interactive reporting for descriptor_exploration; no statistics are recomputed here."""
import base64
import html
import json


SCRIPT = r'''
const data=JSON.parse(document.getElementById('report-data').textContent);
let page=0;
const features=data.features;
const el=(tag,text,cls)=>{const n=document.createElement(tag);if(text!==undefined)n.textContent=text;if(cls)n.className=cls;return n;};
const fmt=x=>x===null||x===undefined?'—':Number(x).toPrecision(4).replace(/\.0+(?=$|e)/,'');
const field=(r,k)=>k==='pearson'||k==='spearman'?Math.abs(r[k]):k==='forest'?r.forest_importance?.mean??-Infinity:r[k];
function selected(){const q=document.getElementById('search').value.toLowerCase(), evidence=document.getElementById('evidence').value;
const sort=document.getElementById('sort').value;
return features.filter(r=>r.name.toLowerCase().includes(q)&&(!evidence||r.evidence===evidence))
.sort((a,b)=>(field(b,sort)-field(a,sort))||a.name.localeCompare(b.name));}
function download(name,rows){const quote=x=>'"'+String(x).replaceAll('"','""')+'"';
const blob=new Blob([rows.map(row=>row.map(quote).join(',')).join('\n')+'\n'],{type:'text/csv;charset=utf-8'});
const a=el('a');a.href=URL.createObjectURL(blob);a.download=name;a.click();setTimeout(()=>URL.revokeObjectURL(a.href),1000);}
function examples(indices){const obs=data.observations;if(!indices.length)return [];
const nat=i=>obs.structures?.[i]?.natoms??Infinity;
let rest=indices.slice(), result=[];const values=indices.map(i=>obs.activity[i]);
const lo=values.reduce((a,b)=>Math.min(a,b),Infinity),hi=values.reduce((a,b)=>Math.max(a,b),-Infinity);
for(const t of [lo,hi,(lo+hi)/2]){if(!rest.length)break;
rest.sort((a,b)=>Math.abs(obs.activity[a]-t)-Math.abs(obs.activity[b]-t)||nat(a)-nat(b)||obs.identifiers[a].localeCompare(obs.identifiers[b]));result.push(rest.shift());}
return result.sort((a,b)=>obs.activity[a]-obs.activity[b]);}
function inspect(indices,caption){const box=document.getElementById('inspect');box.replaceChildren(el('h2',caption),el('p',indices.length+' molecules. Examples span low, middle and high activity; ties favor smaller structures when supplied.'));
const button=el('button','Export all identifiers and activities');button.onclick=()=>download('bucket.csv',[['ID','activity'],...indices.map(i=>[data.observations.identifiers[i],data.observations.activity[i]])]);box.append(button);
const cards=el('div',undefined,'cards');for(const i of examples(indices)){const card=el('div',undefined,'card'),obs=data.observations;
card.append(el('strong',obs.identifiers[i]),el('p','Activity '+fmt(obs.activity[i])));
if(data.depictions?.[i]){const image=el('img');image.src=data.depictions[i];image.alt='Structure of '+obs.identifiers[i];card.append(image);}
else if(obs.structures?.[i])card.append(el('code',obs.structures[i].smiles));
cards.append(card);}box.append(cards);box.scrollIntoView({behavior:'smooth',block:'nearest'});}
function plot(canvas,row){const ctx=canvas.getContext('2d'),w=canvas.width,h=canvas.height;
const points=row.points,curve=row.curve;let xmin=curve[0][0],xmax=curve.at(-1)[0];
let ys=points.map(p=>p[1]).concat(curve.map(p=>p[1]));let ymin=Math.min(...ys),ymax=Math.max(...ys);if(ymax===ymin){ymax+=.5;ymin-=.5;}
const left=65,right=w-20,top=15,bottom=h-45;
const sx=x=>left+(x-xmin)/(xmax-xmin||1)*(right-left),sy=y=>bottom-(y-ymin)/(ymax-ymin)*(bottom-top);
ctx.clearRect(0,0,w,h);ctx.strokeStyle='#8b9299';ctx.beginPath();ctx.moveTo(left,top);ctx.lineTo(left,bottom);ctx.lineTo(right,bottom);ctx.stroke();
ctx.fillStyle='#273444';ctx.font='12px sans-serif';ctx.fillText(fmt(xmin),left,bottom+17);ctx.fillText(fmt(xmax),right-45,bottom+17);ctx.fillText(row.name,w/2,bottom+35);
ctx.fillText(fmt(ymax),3,top+8);ctx.fillText(fmt(ymin),3,bottom);ctx.fillText('activity',3,h/2-15);
ctx.fillStyle='#647e9e55';for(const [x,y] of points){ctx.beginPath();ctx.arc(sx(x),sy(y),2.5,0,7);ctx.fill();}
ctx.strokeStyle='#b43e25';ctx.lineWidth=2;ctx.beginPath();curve.forEach(([x,y],i)=>i?ctx.lineTo(sx(x),sy(y)):ctx.moveTo(sx(x),sy(y)));ctx.stroke();
ctx.fillStyle='#126b48';ctx.strokeStyle='#126b48';for(const b of row.buckets.at(-1).groups){const x=sx((b.low+b.high)/2),y=sy(b.mean);if(b.mean_ci){ctx.beginPath();ctx.moveTo(x,sy(b.mean_ci[0]));ctx.lineTo(x,sy(b.mean_ci[1]));ctx.stroke();}ctx.beginPath();ctx.arc(x,y,4,0,7);ctx.fill();}
canvas.onmousemove=e=>{const rect=canvas.getBoundingClientRect(),x=(e.clientX-rect.left)*w/rect.width,y=(e.clientY-rect.top)*h/rect.height;
let best=null,d=100;for(const p of points){const v=(sx(p[0])-x)**2+(sy(p[1])-y)**2;if(v<d){d=v;best=p;}}
canvas.title=best?data.observations.identifiers[best[2]]+': '+row.name+'='+fmt(best[0])+', activity='+fmt(best[1]):'';};}
function details(row){const detail=el('details'),summary=el('summary','Curve, buckets, correlated properties and molecule exports');detail.append(summary);
if(row.notes.length)detail.append(el('p',row.notes.join('; '),'warning'));
const canvas=el('canvas');canvas.width=650;canvas.height=285;canvas.setAttribute('aria-label',row.name+' versus measured activity');detail.append(canvas);
const describe=el('p','Gray: observations (possibly reduced for display). Red: descriptive full-data fit. Green: finest retained bucket means with approximate 95% mean intervals. The red curve is not a held-out prediction.');detail.append(describe);
const cv=el('p','Spline gain by repeat: '+row.spline_cv_repeats.map(fmt).join(', ')+'. Linear CV gain: '+fmt(row.linear_cv_gain)+'. IQR fitted activity difference: '+fmt(row.effect_iqr)+'.');detail.append(cv);
for(const profile of row.buckets){detail.append(el('h3',profile.requested+' requested buckets → '+profile.actual+' retained; η² '+fmt(profile.eta_squared)+'; CV gain '+fmt(profile.cv_gain)));
const table=el('table'),head=el('tr');for(const title of ['Descriptor range','n','Mean',...(data.show_median?['Median']:[]),'SD','95% mean interval','Molecules'])head.append(el('th',title));table.append(head);
for(const group of profile.groups){const tr=el('tr');for(const value of [fmt(group.low)+' … '+fmt(group.high),group.n,fmt(group.mean),...(data.show_median?[fmt(group.median)]:[]),fmt(group.sd),group.mean_ci?.map(fmt).join(' … ')??'—'])tr.append(el('td',value));
const cell=el('td'),button=el('button','Inspect / export');button.onclick=()=>inspect(group.members,row.name+': '+fmt(group.low)+' … '+fmt(group.high));cell.append(button);tr.append(cell);table.append(tr);}detail.append(table);}
const near=el('p','Correlated descriptors (pairwise complete Spearman): '+(row.neighbors.map(n=>n.name+' ('+fmt(n.spearman)+', n='+n.n+')').join('; ')||'None above the configured threshold.'));detail.append(near);
detail.addEventListener('toggle',()=>{if(detail.open)plot(canvas,row);});return detail;}
function render(){const rows=selected(),size=data.page_size;page=Math.max(0,Math.min(page,Math.ceil(rows.length/size)-1));
const body=document.getElementById('rows');body.replaceChildren();for(const row of rows.slice(page*size,(page+1)*size)){
const tr=el('tr');tr.append(el('td',row.name),el('td',fmt(row.score)),el('td',row.evidence,row.evidence));
for(const value of [row.n,row.distinct,row.pearson,row.spearman,row.mutual_information,row.spline_cv_gain,row.curve_range,row.familywise_p,row.forest_importance?.mean])tr.append(el('td',fmt(value)));
body.append(tr);const lower=el('tr'),cell=el('td');cell.colSpan=12;cell.append(details(row));lower.append(cell);body.append(lower);}
document.getElementById('count').textContent=rows.length+' of '+features.length+' retained descriptors; page '+(page+1)+' of '+Math.max(1,Math.ceil(rows.length/size));
document.getElementById('prev').disabled=page===0;document.getElementById('next').disabled=(page+1)*size>=rows.length;}
for(const id of ['search','sort','evidence'])document.getElementById(id).addEventListener('input',()=>{page=0;render();});
document.getElementById('prev').onclick=()=>{page--;render();};document.getElementById('next').onclick=()=>{page++;render();};
document.getElementById('export').onclick=()=>download('descriptor_scores.csv',[
['descriptor','score','evidence','n','Pearson','Spearman','MI','spline_CV_gain','curve_range','familywise_p','forest_importance'],
...selected().map(r=>[r.name,r.score,r.evidence,r.n,r.pearson,r.spearman,r.mutual_information,r.spline_cv_gain,r.curve_range,r.familywise_p??'',r.forest_importance?.mean??''])]);
render();
'''


def representative_indices(groups, observations):
    """Select image candidates using the same endpoint/midpoint rule as the browser."""
    result = []
    activity = observations['activity']
    structures = observations.get('structures', [])
    for group in groups:
        remaining = list(group['members'])
        if not remaining:
            continue
        low = min(activity[i] for i in remaining)
        high = max(activity[i] for i in remaining)
        for target in (low, high, low * 0.5 + high * 0.5):
            if not remaining:
                break
            def key(i):
                structure = structures[i] if structures else None
                return (abs(activity[i] - target), structure['natoms'] if structure else float('inf'),
                        observations['identifiers'][i])
            chosen = min(remaining, key=key)
            remaining.remove(chosen)
            result.append(chosen)
    return result


def build_html(data, args):
    rows = data['features'][:args.top_n] if args.top_n else data['features']
    payload = dict(data, features=rows, page_size=args.page_size, show_median=not args.nomedian)
    # Embed only bounded example images. Statistical results and molecule exports
    # remain complete; SMILES are shown when an image was not embedded.
    structures = data['observations'].get('structures')
    images = {}
    if structures:
        from activity_report import depiction_svg
        candidates = representative_indices([g for row in rows for profile in row['buckets']
                                              for g in profile['groups']], data['observations'])
        attempted = set()
        for i in candidates:
            if i in attempted or not structures[i]:
                continue
            if len(attempted) >= args.max_structures:
                break
            attempted.add(i)
            svg = depiction_svg(structures[i]['smiles'])
            if svg:
                images[str(i)] = 'data:image/svg+xml;base64,' + base64.b64encode(svg.encode()).decode()
    payload['depictions'] = images
    encoded = json.dumps(payload, allow_nan=False).replace('<', '\\u003c').replace('&', '\\u0026')
    metadata = data['metadata']
    warnings = ''.join('<li>' + html.escape(text) + '</li>' for text in metadata['warnings'])
    omitted = ''.join('<tr><td>' + html.escape(row['name']) + '</td><td>' + html.escape(row['reason']) +
                      '</td></tr>' for row in data['omitted'])
    forest = data.get('forest')
    forest_note = f"Forest held-out gain: {forest['cv_gain']:.3f}." if forest else 'Forest analysis disabled.'
    description = (f"{metadata['matched']} molecules; {metadata['analysed_descriptors']} descriptors analysed; "
                   f"{len(rows)} embedded in this HTML. Activity: {metadata['activity_column']}; "
                   f"transform {metadata['transform']}; desirable direction {metadata['direction']}; units {metadata['units'] or 'unspecified'}.")
    return '''<!doctype html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1"><title>''' + html.escape(args.title) + '''</title>
<style>
body{font:15px system-ui,sans-serif;color:#233144;background:#f5f7fa;margin:25px}h1{font-size:26px}
section,details{background:white;padding:12px;margin:10px 0;border:1px solid #dde3ea;border-radius:6px}
table{border-collapse:collapse;width:100%;background:white}th,td{padding:8px;border-bottom:1px solid #dde3ea;text-align:left;vertical-align:top}th{background:#e9eef5;white-space:nowrap}
.controls{display:flex;gap:12px;flex-wrap:wrap;align-items:center;margin:15px 0}input,select,button{font:inherit;padding:6px}button,summary{cursor:pointer}
.warning,.weak{color:#95451a}.tentative{color:#876000}.supported{color:#17693c}canvas{max-width:100%;background:white}
.cards{display:flex;gap:15px;flex-wrap:wrap}.card{padding:10px;max-width:250px;border:1px solid #ddd}.card img{display:block;width:220px;height:160px}.card code{overflow-wrap:anywhere}
.scroll{overflow-x:auto}code{font-size:12px}small{color:#536170}
</style></head><body><h1>''' + html.escape(args.title) + '''</h1><p>''' + html.escape(description) + '''</p>
<ul class="warning">''' + warnings + '''</ul><p>''' + html.escape(forest_note) + '''</p>
<details><summary>How to read this report</summary>
<p>Scores prioritize inspection; they are not probabilities, causal effects, or evidence that changing a property will change an outcome. A high rank can still have weak evidence. Positive and negative correlation are equally eligible for a high score.</p>
<p>CV gain is 1 minus held-out squared error divided by the squared error of the training-fold mean. Zero gives no improvement; negative is worse. Reported spline complexity is fixed. Curves use all observed data for description only. Mean intervals are approximate t intervals assuming independent molecules; they are not simultaneous intervals and do not account for chemical-series dependence.</p>
<p>Supported requires a familywise shuffle p ≤0.05, positive spline gain in every repeat, mean gain at least the configured threshold, sufficient descriptive effect, and at least 30 observations. This is exploratory calibration on this dataset, not independent confirmation. With 19 shuffles, the finest possible p is 0.05. Forest importance can be diluted by correlated substitutes.</p>
<p>Composite formula: <code>''' + html.escape(metadata['score_definition']) + '''</code>. The MI transformation only puts it on a bounded heuristic scale. Bucket η² and MI enter one dependence component. Raw metrics remain available. Sparse adjacent buckets are merged and ties are never split. JSON contains all analysed descriptors; the HTML limit applies before interactive filtering.</p>
</details>
<div class="controls"><label>Filter <input id="search" type="search"></label>
<label>Sort <select id="sort"><option value="score">Composite</option><option value="spline_cv_gain">Spline CV gain</option><option value="pearson">|Pearson|</option><option value="spearman">|Spearman|</option><option value="mutual_information">Mutual information</option><option value="curve_range">Descriptive effect range</option><option value="forest">Forest importance</option></select></label>
<label>Evidence <select id="evidence"><option value="">All</option><option>supported</option><option>tentative</option><option>weak</option></select></label>
<button id="export">Export displayed ranking</button></div>
<p id="count"></p><button id="prev">Previous</button> <button id="next">Next</button>
<div class="scroll"><table><thead><tr><th>Descriptor</th><th>Score</th><th>Evidence</th><th>n</th><th>Distinct</th><th>Pearson</th><th>Spearman</th><th>MI</th><th>Spline CV gain</th><th>Curve range</th><th>Familywise p</th><th>Forest importance</th></tr></thead><tbody id="rows"></tbody></table></div>
<section id="inspect"><p>Open a descriptor and select a bucket to inspect examples or export its molecule identifiers.</p></section>
<details><summary>Omitted descriptors (''' + str(len(data['omitted'])) + ''')</summary><table>''' + omitted + '''</table></details>
<script id="report-data" type="application/json">''' + encoded + '''</script><script>''' + SCRIPT + '''</script></body></html>'''
