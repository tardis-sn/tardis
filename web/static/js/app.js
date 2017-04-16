var jsonConf = {};
    
var downloadConfFile = function(jsonData){    
    yaml_data = YAML.stringify(jsonData,5);
    // console.log(jsonConf)
    // console.log(yaml_data)
    MIME_TYPE = "text/yml";
    theBlob = new Blob([yaml_data], {type: MIME_TYPE});
    a = document.getElementById("dwnldYamlFile");
    a.download = "tardis_config.yml";
    a.href = window.URL.createObjectURL(theBlob);
    a.dataset.downloadurl = [MIME_TYPE, a.download, a.href].join(":");
};

var changeType = function(val,type){
    switch (type) {
        case "int":
            return parseInt(val)
            break; 
        case "float":
            return parseFloat(val)
            break; 
        default: 
            return val
    }
}

var manipulateValue = function(field,parentTagsStr,value,valType,remFlag){

    if (parentTagsStr === ""){
        temp = jsonConf;
    }else{    
        parentTags = parentTagsStr.split('-');
        depth = parentTags.length;                
        if(!jsonConf.hasOwnProperty(parentTags[0])){
            jsonConf[parentTags[0]] = {}
        }
        temp = jsonConf[parentTags[0]]
        for(var x=1; x <= parentTags.length-1 ; x++){
            if(!temp.hasOwnProperty(parentTags[x])){
                temp[parentTags[x]] = {}
            }
            temp = temp[parentTags[x]]
        }
    }

    if(remFlag){
        delete temp[field]
    }else{
        if(valType === "int" || valType === "float"){
            value = changeType(value,valType)
        }
        temp[field] = value;
    }

    console.log(jsonConf)
}

$(document).ready(function(){
    $('.badge').tooltip();
    
    $(".containerProp").on('change', function() {
        $($("option:not(:selected)", this)).each(function() {
            $($(this).attr("data-targetCont")).slideUp();
        });
        $($("option:selected", this).attr("data-targetCont")).slideDown();
    });
    
    $('input[type=text]').not('.tags-inp').on('change',function(){
        if($(this).val() === ""){
            remFlag = true;
        }else{
            remFlag = false;
        }
        manipulateValue($(this).attr('data-field'), $(this).attr('data-parentTags'), $(this).val(), $(this).attr('data-type'), remFlag)
    })

    $('input.tags-inp').on('itemAdded', function(event) {
        tag = event.item;
        manipulateValue(tag.split(':')[0], $(this).attr('data-parentTags'), tag.split(':')[1], $(this).attr('data-type'), false)
    });

    $('input.tags-inp').on('itemRemoved', function(event) {
        tag = event.item;
        manipulateValue(tag.split(':')[0], $(this).attr('data-parentTags'), tag.split(':')[1], $(this).attr('data-type'), true)
    });

    $('input[type=checkbox]').on('change',function(){
        manipulateValue($(this).attr('data-field'), $(this).attr('data-parentTags'), $(this).prop("checked"), $(this).attr('data-type'), false)
    })

    $('select').not('.containerProp').on('change',function(){
        if($(this).val() === ""){
            remFlag = true;
        }else{
            remFlag = false;
        }
        manipulateValue($(this).attr('data-field'), $(this).attr('data-parentTags'), $(this).val(), $(this).attr('data-type'), remFlag)
    });

    $('select.containerProp').on('change',function(){
        parenTags = $(this).attr('data-parentTags').split('-')
        manipulateValue(parenTags.pop(), parenTags.join('-'), "", $(this).attr('data-type'), true);
        if($(this).val() !== ""){
            manipulateValue($(this).attr('data-field'), $(this).attr('data-parentTags'), $(this).val(), $(this).attr('data-type'), false)
        }
        
    });
    
    $('.nav li a').on('click',function(){
        if($(this).parent().is(':not(.active)')){
            activeDiv = $('.nav li.active a').attr('data-target');
            $('.nav li').removeClass('active');
            targetDiv = $(this).attr('data-target');
            $(this).parent().addClass('active');
            $(activeDiv).fadeOut('medium',function(){
                $(targetDiv).fadeIn();
            })
        }
    });
    
    $('#dwnldYamlFile').on('click',function(){
        downloadConfFile(jsonConf);
    });
    
})