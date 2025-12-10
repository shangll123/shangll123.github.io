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
			Shang Lab
		</div>
	</div>	
</div> 




<img src="/assets/themes/lab/images/logo/lulu.jpeg" alt="photo" width="300" class="center">

<br clear="left"/>
<hr/>

[[Publications]](/papers/) 
[[Google Scholar<i class="ai ai-google-scholar"></i>]](https://scholar.google.com/citations?user=7FEgLPkAAAAJ&hl=en&authuser=1) [[GitHub<i class="fa fa-github"></i>]](https://github.com/shangll123) [[Twitter<i class="fa fa-twitter"></i>]](https://twitter.com/shang_lulu).


<hr/>

**About Me**:

My name is Lulu Shang. I am a tenure-track assistant professor in the Department of Biostatistics at MD Anderson Cancer Center, and started Nov 2023. Previously, I received my Ph.D. from the Department of Biostatistics at the University of Michigan, under supervision of Prof. [Xiang Zhou](http://xzlab.org). Prior to that, I obtained my bachelor’s degree in Biology from the Zhiyuan Honored College at Shanghai Jiao Tong University.
 
<hr/>

**Research Interests**:

My research program centers on developing innovative statistical and computational methods for spatial biology, single-cell genomics, and integrative multi-omics. I design scalable, interpretable, and data-driven models that capture complex molecular, cellular, and spatial organization in human tissues, with the broader goal of advancing our understanding of disease mechanisms across diverse biological systems. By integrating statistical modeling, deep learning, and high-dimensional genomic profiling, my lab builds tools that enable quantitative characterization of tissue architecture, cellular heterogeneity, and gene regulatory processes. We work closely with experimental and clinical collaborators to ensure that our methods generate biologically meaningful and clinically actionable insights, including applications in cancer when improved spatial and molecular resolution can inform diagnosis, prognosis, and therapeutic strategies. Ultimately, my research aims to create broadly applicable computational frameworks that strengthen the bridge between modern genomics technologies and biomedical discovery.

**Keywords**:

- *Statistical Methods*: High-dimensional data analysis, mixed effects models, longitudinal data analysis, spatial statistics, multiple instance learning, 
non-parametric models, dimension reduction, data integration, statistical computing, network analysis, and deep learning.
- *Applications*: High dimensional genetic and genomic data including spatial transcriptomics, single cell RNA sequencing studies, spatial and single cell multi-omics, digital pathology, genome-wide association studies, genome-wide quantitative loci mapping, methylation studies, and gene regulatory network.


<hr/>

**Funding**
- Institutional Research Grant, $75k direct cost, 2/1/2025-1/31/2026 (Awarded to top 20 proposals among all faculty applicants at MD Anderson)

<hr/>

**Selected Awards**
- IMS New Researcher Travel Grant Award, Institute of Mathematical Statistics, 2024
- ProQuest Distinguished Dissertation Awards, Rackham Graduate School, University of Michigan, 2024 (Awarded to top 10 Ph.D. dissertations among 800+ applicants, the first recipient in the Department of Biostatistics since 2007)
- Charles J. Epstein Trainee Awards, Predoctoral semifinalist, American Society of Human Genetics, 2022
- Reviewer’s Choice Poster Award, American Society of Human Genetics, 2021
- Excellence in Research Award Honorable Mention, Department of Biostatistics, University of Michigan, 2020
- Most Interesting Methodological Advancement Poster Award, MIDAS Annual Symposium, University of Michigan, 2019

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



