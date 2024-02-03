---
layout: default
title: Lulu Shang
categories:
 - home
---
{% include JB/setup %}
{% for page in site.categories.misc %}
{% if page.homepage %}
	{% assign homepage = page %}
{% endif %}
{% endfor %}

<link rel="stylesheet" href="https://cdn.jsdelivr.net/gh/jpswalsh/academicons@1/css/academicons.min.css">

<div class="row">
	<div class="col-md-12">
		<!-- <object class="pull-left biglogo" data="assets/themes/lab/images/logo/logo-none.svg" type="image/svg+xml"></object> -->
			<div class="bigtitle logobox">	
			Lulu Shang
		</div>
	</div>	
</div> 





<br clear="left"/>
<hr/>

[[Publications]](/papers/) 
[[Google Scholar<i class="ai ai-google-scholar"></i>]](https://scholar.google.com/citations?hl=en&user=tkt5ZOYAAAAJ&view_op=list_works&sortby=pubdate) [[GitHub<i class="fa fa-github"></i>]](https://github.com/shangll123) [[Twitter<i class="fa fa-twitter"></i>]](https://twitter.com/shang_lulu).


<hr/>
I am a tenure-track assistant professor in the Department of Biostatistics at MD Anderson Cancer Center. I received my PhD degree in the Department of Biostatistics at the University of Michigan. I have been working with Prof. [Xiang Zhou](http://xzlab.org) on developing statistical methods for genetic and genomic datasets. My recent focus is on developing methods for analyzing spatial transcriptomics. I also work closely with Prof. [Jennifer Smith](https://sph.umich.edu/faculty-profiles/smith-jennifer.html) in the Department of Epidemiology at the University of Michigan on large-scale quantitative trait mapping in African Americans in the GENOA study.
	


	
<hr/>


**Research Interests**:

My research interest is in developing effective and efficient statistical and computational methods in single cell and spatial transcriptomics to address critical biological problems such as omics data dimension reduction and integration. My specific focus includes: 1). developing in statistical and computational methods for data analysis in single-cell RNA-seq and spatial transcriptomics; 2). Bridging between population-level genome-wide association studies (GWASs) and individual-level single-cell and spatial transcriptomics data to unravel the underlying mechanisms of disease etiology; 3). Integrating multi-omics data with patient outcomes to foster translational research, ultimately connecting molecular insights with clinical applications. 

**Keywords**:

- *Statistical Methods*: High-dimensional data analysis, mixed effects models, longitudinal data analysis, spatial statistics,
non-parametric models, dimension reduction, data integration, statistical computing, network analysis and machine learning.
- *Applications*: High dimensional genetic and genomic data including spatial transcriptomics, single cell RNA sequencing studies, spatial and single cell multi-omics, genome-wide association studies, genome-wide quantitative
loci mapping, methylation studies, and gene regulatory network.

**Open positions**:

Applications are invited for a postdoctoral fellow position in my research group. The successful candidates will be working on various research topics in developing statistical methods and computational tools in the field of single cell and spatial transcriptomics. The successful candidates will be offered with competitive benefits and have the opportunity to analyze a variety of large-scale data types. Applicants should have, or be studying for, a PhD in biostatistics, statistics, computer science, bioinformatics, computational biology, mathematics, or related quantitative discipline. A strong computational background is preferred. Applicants should send a CV, a short statement of research interests, and contact information of three referees to: Lulu Shang [lshang@mdanderson.org](lshang@mdanderson.org). Review of applications will begin immediately and continue until the position is filled.

We also have open positions for Research Assistant. The Research Assistant positions can be offered to students in MD Anderson, Rice, and UTHealth. Please feel free to reach out [lshang@mdanderson.org](lshang@mdanderson.org) if you are interested in joining our lab.

<hr/>



<br />

<hr/>

<div class="row">
	<div class="col-md-12">
		<div class="head">
			{{ homepage.content }}
		</div>
	</div>				
</div>

<div class="row">
	

	
	<div class="col-md-4">
		<div class="head">
			<a class="off" href="/papers/">Recent Papers
			</a>
		</div>
		<div class="bigspacer"></div>
		<div class="feedbox pad-left">		
			{% for paper in site.categories.papers limit:10 %}
				<div class="note-title">
					<i class="fa fa-file-text-o fa-fw"></i>
					<a class="off" href="{{ paper.url }}">
					{{ paper.title }}
					</a>
					<br/>
					<div class='shortref note'>
					{{ paper.shortref }}
					</div>
				</div>
				<div class="smallspacer"></div>
				<div class="smallnote">
					Published
					{{ paper.date | date_to_string }}
				</div>
				<div class="spacer"></div>	
				<div class="spacer"></div>				
			{% endfor %}
		</div>
		<div class="bigspacer"></div>		
	</div>
	
    	<div class="col-md-4">
		<div class="head">
			<a class="off" href="/news/">News</a>
		</div>
		<div class="bigspacer"></div>
		<div class="feedbox pad-left">
			{% for news in site.categories.news limit:5 %}
			
				{% for member in site.categories.team %}
					{% if member.handle == news.author_handle %}
						{% assign author = member %}
					{% endif %}
				{% endfor %}		
				
				<div class="note-title">
					<i class="fa fa-bullhorn"></i>
					<a class="off" href="{{ news.url }}">
					{{ news.title }}
					</a>
				</div>
				<div class="note">
					{{ news.content }}
				</div>
				<div class="smallspacer"></div>
				<div class="smallnote">
					Posted
					{{ news.date | date_to_string }}
					{% if author %}
					by <a class="off" href="{{ author.url }}">{{ author.handle }}</a>
					{% endif %}						
				</div>
				<div class="spacer"></div>	
				<div class="spacer"></div>				
			{% endfor %}
		</div>
		<div class="bigspacer"></div>		
	</div>

	<div class="col-md-4">
		<div class="head">
			<a class="off" href="/blog/">Blogs</a>
		</div>
		<div class="bigspacer"></div>
		<div class="feedbox pad-left">
			{% for blog in site.categories.blog limit:5 %}
				<div class="note-title">
					<i class="fa fa-comment-o fa-fw"></i>
					<a class="off" href="{{ blog.url }}">
					{{ blog.title }}
					</a>
					<br/>
					<div class='shortref note'>
					{{ blog.description }}
					</div>
				</div>
				
				<div class="smallspacer"></div>
				<div class="smallnote">
					Posted
					{{ blog.date | date_to_string }}
					{% if author %}
					by <a class="off" href="{{ author.url }}">{{ author.handle }}</a>
					{% endif %}						
				</div>
				<div class="spacer"></div>	
				<div class="spacer"></div>
				
				<div class="smallspacer"></div>
				<div class="spacer"></div>
				<div class="spacer"></div>
			{% endfor %}
		</div>
		<div class="bigspacer"></div>
	</div>
	
	<!-- <div class="col-md-4">
		<div class="head">
			<a class="off" href="/projects/">Projects</a>
		</div>
		<div class="bigspacer"></div>
		<div class="feedbox pad-left">
			{% for project in site.categories.project limit:4 %}
				<div class="note-title">
					<i class="fa fa-terminal"></i>
					<a class="off" href="{{ project.url }}">
					{{ project.title }}
					</a>
					<br/>
					<div class='shortref note'>
					{{ project.tags }}
					</div>
				</div>
				<div class="smallspacer"></div>
				<div class="spacer"></div>
				<div class="spacer"></div>
			{% endfor %}
		</div>
		<div class="bigspacer"></div>
	</div> -->


</div>

<div class="bigspacer"></div>

<img src="/assets/themes/lab/images/logo/background2.png" alt="photo" width="800">



